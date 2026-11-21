# Made by Debarshi Bose (with help from Gemini)

import pandas as pd
import matplotlib.pyplot as plt

df = pd.read_csv('analysis.csv')

fig, ax = plt.subplots(figsize=(12, 8))

sc = ax.scatter(df['Q_Fraction'], df['Energy'],
                c=df['Step'], cmap='viridis_r',
                s=25, alpha=0.7, edgecolors='none')

ax.scatter(df['Q_Fraction'].iloc[0], df['Energy'].iloc[0],
           color='white', s=200, edgecolor='black', linewidth=1.5,
           label='Unfolded State (Start)', zorder=5)

ax.scatter(df['Q_Fraction'].iloc[-1], df['Energy'].iloc[-1],
           color='black', s=200, edgecolor='white', linewidth=1.5,
           label='Native State\n(79% of actual native crystallographic conformation)', zorder=5)

ax.set_xlabel('Fraction of Native Contacts (Q)', fontsize=14, fontweight='bold')
ax.set_ylabel('System Energy', fontsize=14, fontweight='bold')
ax.set_title('2D Conformational Free Energy Landscape', fontsize=18, fontweight='bold', pad=15)

ax.set_xlim(-0.05, 1.05)

cbar = plt.colorbar(sc, ax=ax, pad=0.02)
cbar.set_label('Monte Carlo Step (Time)', rotation=270, labelpad=20, fontsize=12, fontweight='bold')

ax.grid(True, linestyle='--', alpha=0.5)
ax.legend(fontsize=11, loc='upper right')

ax.annotate('Hydrophobic Collapse\n(Rapid energy drop)',
            xy=(0.15, df['Energy'].max() - 2500),
            xytext=(0.20, df['Energy'].max() - 1000),
            arrowprops=dict(facecolor='black', arrowstyle='->', lw=1.5),
            fontsize=12, bbox=dict(facecolor='white', alpha=0.8, edgecolor='none'))

plt.savefig('2d_landscape.png', dpi=300, bbox_inches='tight')
print("Rectified 2D landscape generated successfully.")