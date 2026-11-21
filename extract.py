# Made by Debarshi Bose (with help from Gemini)

import math
import csv


def get_dihedral(p0, p1, p2, p3):
    b1 = [p1[0]-p0[0], p1[1]-p0[1], p1[2]-p0[2]]
    b2 = [p2[0]-p1[0], p2[1]-p1[1], p2[2]-p1[2]]
    b3 = [p3[0]-p2[0], p3[1]-p2[1], p3[2]-p2[2]]

    def cross(a, b):
        return [a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2], a[0]*b[1]-a[1]*b[0]]

    def dot(a, b):
        return sum(x*y for x, y in zip(a, b))

    def mag(a):
        return math.sqrt(dot(a, a))

    n1, n2 = cross(b1, b2), cross(b2, b3)
    m1 = cross(n1, b2)
    return math.atan2(dot(m1, n2) / mag(b2), dot(n1, n2))


energies = {}
with open('analysis.csv', 'r') as f:
    reader = csv.reader(f)
    next(reader)
    for row in reader:
        energies[int(row[0])] = float(row[1])

with open('trajectory.pdb', 'r') as f_in, open('torsion_data.csv', 'w') as f_out:
    f_out.write('Step,Energy,Torsion1,Torsion2\n')
    step = 0
    coords = []

    for line in f_in:
        if line.startswith('REMARK MC Step'):
            step = int(line.split()[-1])
            coords = []
        elif line.startswith('ATOM'):
            x, y, z = float(line[30:38]), float(line[38:46]), float(line[46:54])
            coords.append((x, y, z))
        elif line.startswith('ENDMDL'):
            if step in energies and len(coords) > 5:
                mid = len(coords) // 2
                t1 = get_dihedral(coords[mid-2], coords[mid-1], coords[mid], coords[mid+1])
                t2 = get_dihedral(coords[mid-1], coords[mid], coords[mid+1], coords[mid+2])
                f_out.write(f"{step},{energies[step]},{t1},{t2}\n")

print("Extraction complete: torsion_data.csv generated.")