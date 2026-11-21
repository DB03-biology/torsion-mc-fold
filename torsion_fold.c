/* Made by Debarshi Bose (with help from Gemini) */

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <ctype.h>
#include <time.h>

#define SAVE_FREQ 500000
#define TEMP 1.35
#define K_BOND 100.0
#define K_ANGLE 20.0
#define K_DIH 15.0
#define EPSILON_NAT 1.5
#define CONTACT_CUTOFF 8.0
#define MAX_DISP 0.10
#define TARGET_Q 0.95
#define MAX_STEPS 1500000000

typedef struct { double x, y, z; } Vector3D;

typedef struct { char sym; char name[4]; double radius; double hydrophobicity; } AAParam;

AAParam aa_db[20] = {
    {'A', "ALA", 5.0,  1.8}, {'R', "ARG", 6.5, -4.5}, {'N', "ASN", 5.4, -3.5},
    {'D', "ASP", 5.4, -3.5}, {'C', "CYS", 5.4,  2.5}, {'Q', "GLN", 5.9, -3.5},
    {'E', "GLU", 5.9, -3.5}, {'G', "GLY", 4.3, -0.4}, {'H', "HIS", 6.0, -3.2},
    {'I', "ILE", 6.0,  4.5}, {'L', "LEU", 6.0,  3.8}, {'K', "LYS", 6.4, -3.9},
    {'M', "MET", 6.2,  1.9}, {'F', "PHE", 6.4,  2.8}, {'P', "PRO", 5.5, -1.6},
    {'S', "SER", 5.1, -0.8}, {'T', "THR", 5.6, -0.7}, {'W', "TRP", 6.8, -0.9},
    {'Y', "TYR", 6.5, -1.3}, {'V', "VAL", 5.9,  4.2}
};

double *r0_bond;
double *th0_angle;
double *dih0_angle;
double **r0_cmap;
int **cmap;
int total_native_contacts = 0;

double calc_local_energy(Vector3D coords[], int seq_idx[], int atom, int n_res);
double calc_total_energy(Vector3D coords[], int seq_idx[], int n_res);
double calc_Q(Vector3D curr[], int n_res);
int parse_native_pdb(const char *filename, Vector3D coords[], int max_res);
void write_pdb_frame(FILE *fp, Vector3D coords[], int seq_idx[], int n_res, int step, int frame);

int get_aa_idx(char sym) {
    for (int i = 0; i < 20; i++) {
        if (aa_db[i].sym == sym) return i;
    }
    return 0;
}

double dist(Vector3D a, Vector3D b) {
    return sqrt((a.x-b.x)*(a.x-b.x) + (a.y-b.y)*(a.y-b.y) + (a.z-b.z)*(a.z-b.z));
}

double dot(Vector3D v1, Vector3D v2) {
    return v1.x*v2.x + v1.y*v2.y + v1.z*v2.z;
}

double mag(Vector3D v) {
    return sqrt(v.x*v.x + v.y*v.y + v.z*v.z);
}

Vector3D cross(Vector3D v1, Vector3D v2) {
    return (Vector3D){v1.y*v2.z - v1.z*v2.y, v1.z*v2.x - v1.x*v2.z, v1.x*v2.y - v1.y*v2.x};
}

double angle(Vector3D a, Vector3D b, Vector3D c) {
    Vector3D v1 = {a.x - b.x, a.y - b.y, a.z - b.z};
    Vector3D v2 = {c.x - b.x, c.y - b.y, c.z - b.z};
    double cos_th = dot(v1, v2) / (mag(v1) * mag(v2));
    if (cos_th > 1.0) cos_th = 1.0;
    if (cos_th < -1.0) cos_th = -1.0;
    return acos(cos_th);
}

double dihedral(Vector3D p0, Vector3D p1, Vector3D p2, Vector3D p3) {
    Vector3D b1 = {p1.x-p0.x, p1.y-p0.y, p1.z-p0.z};
    Vector3D b2 = {p2.x-p1.x, p2.y-p1.y, p2.z-p1.z};
    Vector3D b3 = {p3.x-p2.x, p3.y-p2.y, p3.z-p2.z};

    Vector3D n1 = cross(b1, b2);
    Vector3D n2 = cross(b2, b3);
    Vector3D m1 = cross(n1, b2);

    double x = dot(n1, n2);
    double y = dot(m1, n2) / mag(b2);
    return atan2(y, x);
}

double rand_double() {
    return 2.0 * ((double)rand() / RAND_MAX) - 1.0;
}

int main(int argc, char *argv[]) {
    if (argc < 3) {
        printf("Usage: ./torsion_fold <SEQUENCE> <target_structure.pdb>\n");
        return 1;
    }

    char *seq = argv[1];
    int n_res = strlen(seq);
    char *pdb_filename = argv[2];

    int *seq_idx = (int *)malloc(n_res * sizeof(int));
    for (int i = 0; i < n_res; i++) seq_idx[i] = get_aa_idx(toupper(seq[i]));

    Vector3D *native_coords = (Vector3D *)malloc(n_res * sizeof(Vector3D));
    Vector3D *current_coords = (Vector3D *)malloc(n_res * sizeof(Vector3D));

    int pdb_res = parse_native_pdb(pdb_filename, native_coords, n_res);
    if (pdb_res != n_res) {
        printf("Error: Sequence length (%d) != PDB CA atoms (%d).\n", n_res, pdb_res);
        return 1;
    }

    r0_bond = (double *)malloc((n_res - 1) * sizeof(double));
    th0_angle = (double *)malloc((n_res - 2) * sizeof(double));
    dih0_angle = (double *)malloc((n_res - 3) * sizeof(double));
    cmap = (int **)malloc(n_res * sizeof(int *));
    r0_cmap = (double **)malloc(n_res * sizeof(double *));

    for (int i = 0; i < n_res; i++) {
        cmap[i] = (int *)calloc(n_res, sizeof(int));
        r0_cmap[i] = (double *)calloc(n_res, sizeof(double));
        current_coords[i].x = i * 3.8;
        current_coords[i].y = 0.0;
        current_coords[i].z = 0.0;
    }

    for (int i = 0; i < n_res - 1; i++) r0_bond[i] = dist(native_coords[i], native_coords[i+1]);
    for (int i = 1; i < n_res - 1; i++) th0_angle[i-1] = angle(native_coords[i-1], native_coords[i], native_coords[i+1]);
    for (int i = 0; i < n_res - 3; i++) dih0_angle[i] = dihedral(native_coords[i], native_coords[i+1], native_coords[i+2], native_coords[i+3]);

    for (int i = 0; i < n_res - 3; i++) {
        for (int j = i + 3; j < n_res; j++) {
            double r = dist(native_coords[i], native_coords[j]);
            if (r < CONTACT_CUTOFF) {
                cmap[i][j] = 1;
                cmap[j][i] = 1;
                r0_cmap[i][j] = r;
                r0_cmap[j][i] = r;
                total_native_contacts++;
            }
        }
    }

    srand((unsigned int)time(NULL));
    FILE *out_pdb = fopen("trajectory.pdb", "w");
    FILE *out_csv = fopen("analysis.csv", "w");
    fprintf(out_csv, "Step,Energy,Q_Fraction\n");

    double total_energy = calc_total_energy(current_coords, seq_idx, n_res);
    int accepted = 0, frame = 1, step = 0;
    double current_q = 0.0;

    printf("Folding %d residues. Target Q: %.2f\n", n_res, TARGET_Q);

    while (current_q < TARGET_Q && step < MAX_STEPS) {
        step++;
        int atom = rand() % n_res;
        Vector3D old_pos = current_coords[atom];

        double e_local_old = calc_local_energy(current_coords, seq_idx, atom, n_res);

        current_coords[atom].x += rand_double() * MAX_DISP;
        current_coords[atom].y += rand_double() * MAX_DISP;
        current_coords[atom].z += rand_double() * MAX_DISP;

        double e_local_new = calc_local_energy(current_coords, seq_idx, atom, n_res);
        double delta_e = e_local_new - e_local_old;

        if (delta_e < 0.0 || (exp(-delta_e / TEMP) > ((double)rand() / RAND_MAX))) {
            total_energy += delta_e;
            accepted++;
        } else {
            current_coords[atom] = old_pos;
        }

        if (step % SAVE_FREQ == 0) {
            current_q = calc_Q(current_coords, n_res);
            fprintf(out_csv, "%d,%.3f,%.4f\n", step, total_energy, current_q);
            write_pdb_frame(out_pdb, current_coords, seq_idx, n_res, step, frame++);

            printf("Step %8d | Acc: %4.1f%% | E: %8.1f | Q: %.3f / %.2f\n",
                   step, 100.0 * accepted / step, total_energy, current_q, TARGET_Q);
        }
    }

    if (current_q >= TARGET_Q) {
        printf("\nTarget Fold Reached! Simulation matched %.1f%% of native topology.\n", current_q * 100.0);
    } else {
        printf("\nHit step limit. Simulation matched %.1f%% of native topology.\n", current_q * 100.0);
    }

    fclose(out_pdb);
    fclose(out_csv);
    return 0;
}

double calc_local_energy(Vector3D coords[], int seq_idx[], int atom, int n_res) {
    double energy = 0.0;

    if (atom > 0) energy += K_BOND * pow(dist(coords[atom-1], coords[atom]) - r0_bond[atom-1], 2);
    if (atom < n_res - 1) energy += K_BOND * pow(dist(coords[atom], coords[atom+1]) - r0_bond[atom], 2);

    if (atom > 1) energy += K_ANGLE * pow(angle(coords[atom-2], coords[atom-1], coords[atom]) - th0_angle[atom-2], 2);
    if (atom > 0 && atom < n_res - 1) energy += K_ANGLE * pow(angle(coords[atom-1], coords[atom], coords[atom+1]) - th0_angle[atom-1], 2);
    if (atom < n_res - 2) energy += K_ANGLE * pow(angle(coords[atom], coords[atom+1], coords[atom+2]) - th0_angle[atom], 2);

    for (int k = atom - 3; k <= atom; k++) {
        if (k >= 0 && k < n_res - 3) {
            double dih = dihedral(coords[k], coords[k+1], coords[k+2], coords[k+3]);
            energy += K_DIH * (1.0 - cos(dih - dih0_angle[k]));
        }
    }

    for (int j = 0; j < n_res; j++) {
        if (abs(atom - j) >= 3) {
            double r = dist(coords[atom], coords[j]);
            if (cmap[atom][j]) {
                double ratio = r0_cmap[atom][j] / r;
                double r6 = pow(ratio, 6);
                energy += EPSILON_NAT * (r6*r6 - 2.0*r6);
            } else {
                double sigma_ij = (aa_db[seq_idx[atom]].radius + aa_db[seq_idx[j]].radius) / 2.0;
                double h_i = aa_db[seq_idx[atom]].hydrophobicity;
                double h_j = aa_db[seq_idx[j]].hydrophobicity;
                double epsilon_ij = 0.3 + 0.05 * (h_i * h_j);
                if (epsilon_ij < 0.1) epsilon_ij = 0.1;
                double r6 = pow(sigma_ij / r, 6);
                energy += 4.0 * epsilon_ij * (r6*r6 - r6);
            }
        }
    }
    return energy;
}

double calc_total_energy(Vector3D coords[], int seq_idx[], int n_res) {
    double energy = 0.0;
    for (int i = 0; i < n_res; i++) energy += calc_local_energy(coords, seq_idx, i, n_res);
    return energy / 2.0;
}

double calc_Q(Vector3D curr[], int n_res) {
    if (total_native_contacts == 0) return 0.0;
    int formed = 0;
    for (int i = 0; i < n_res - 3; i++) {
        for (int j = i + 3; j < n_res; j++) {
            if (cmap[i][j]) {
                double r = dist(curr[i], curr[j]);
                if (r <= r0_cmap[i][j] * 1.2 && r >= r0_cmap[i][j] * 0.8) formed++;
            }
        }
    }
    return (double)formed / total_native_contacts;
}

int parse_native_pdb(const char *filename, Vector3D coords[], int max_res) {
    FILE *fp = fopen(filename, "r");
    if (!fp) {
        printf("Cannot open %s\n", filename);
        exit(1);
    }
    char line[256];
    int count = 0;
    while (fgets(line, sizeof(line), fp)) {
        if (strncmp(line, "ATOM", 4) == 0 && strncmp(line + 13, "CA", 2) == 0) {
            if (count >= max_res) break;
            char x_str[9], y_str[9], z_str[9];
            strncpy(x_str, line + 30, 8); x_str[8] = '\0';
            strncpy(y_str, line + 38, 8); y_str[8] = '\0';
            strncpy(z_str, line + 46, 8); z_str[8] = '\0';
            coords[count].x = atof(x_str);
            coords[count].y = atof(y_str);
            coords[count].z = atof(z_str);
            count++;
        }
    }
    fclose(fp);
    return count;
}

void write_pdb_frame(FILE *fp, Vector3D coords[], int seq_idx[], int n_res, int step, int frame) {
    fprintf(fp, "MODEL %8d\n", frame);
    fprintf(fp, "REMARK MC Step %d\n", step);
    for (int i = 0; i < n_res; i++) {
        fprintf(fp, "ATOM  %5d  CA  %3s A%4d    %8.3f%8.3f%8.3f  1.00  0.00           C  \n",
                i + 1, aa_db[seq_idx[i]].name, i + 1, coords[i].x, coords[i].y, coords[i].z);
    }
    fprintf(fp, "ENDMDL\n");
}