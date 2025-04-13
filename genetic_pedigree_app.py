import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle, Circle, Patch
from google.colab import files

# Upload CSV file
uploaded = files.upload()
filename = list(uploaded.keys())[0]

# Load CSV and clean columns
df = pd.read_csv(filename)
df.columns = df.columns.str.strip()

# Define disorder prediction
def get_disorder_details(row):
    details = []

    cb = row['Genotype_ColorBlindness']
    if 'Y' in cb and 'Xc' in cb:
        details.append("Color Blindness")
    elif cb.count('Xc') == 2:
        details.append("Color Blindness")
    elif 'Xc' in cb:
        details.append("CB Carrier")

    hemo = row['Genotype_Hemophilia']
    if 'Y' in hemo and 'Xh' in hemo:
        details.append("Hemophilia")
    elif hemo.count('Xh') == 2:
        details.append("Hemophilia")
    elif 'Xh' in hemo:
        details.append("Hemophilia Carrier")

    sc = row['Genotype_SickleCell']
    if sc == 'SS':
        details.append("Sickle Cell")
    elif 'S' in sc and 'A' in sc:
        details.append("SC Carrier")

    md = row['Genotype_Myotonic']
    if 'n' in md:
        details.append("Myotonic Dystrophy")

    return details if details else ["Healthy"]

# Assign positions in pedigree
def assign_positions(df):
    pos = {}
    level = 0
    roots = df[pd.isna(df['Parents']) | (df['Parents'] == '')]
    for i, row in enumerate(roots.itertuples()):
        pos[row.Name] = (i * 3, level)

    def place_children(parents, parent_x, parent_level):
        children = df[df['Parents'].fillna('').str.contains(parents)]
        for j, child in enumerate(children.itertuples()):
            x = parent_x + j
            y = parent_level - 1
            pos[child.Name] = (x, y)
            place_children(child.Name, x, y)

    for root in roots.itertuples():
        place_children(root.Name, pos[root.Name][0], level)

    return pos

# Color map per disorder
disorder_colors = {
    "Color Blindness": "orange",
    "CB Carrier": "gold",
    "Hemophilia": "blue",
    "Hemophilia Carrier": "skyblue",
    "Sickle Cell": "red",
    "SC Carrier": "lightcoral",
    "Myotonic Dystrophy": "purple",
    "Healthy": "white"
}

# Draw the pedigree
def draw_pedigree(df):
    pos = assign_positions(df)
    fig, ax = plt.subplots(figsize=(18, 10))

    df['Disorders'] = df.apply(get_disorder_details, axis=1)

    for _, row in df.iterrows():
        x, y = pos[row['Name']]
        disorders = row['Disorders']

        # Determine shape based on sex
        color = disorder_colors.get(disorders[0], "gray")
        if row['Sex'].lower() == 'male':
            shape = Rectangle((x, y), 1, 1, facecolor=color, edgecolor='black', lw=2)
        else:
            shape = Circle((x + 0.5, y + 0.5), 0.5, facecolor=color, edgecolor='black', lw=2)

        ax.add_patch(shape)

        # Text Labels
        label = "\n".join(disorders)
        ax.text(x + 0.5, y - 0.3, row['Name'], ha='center', fontsize=9)
        ax.text(x + 0.5, y - 1.0, label, ha='center', fontsize=7)

        # Draw parent connections
        if pd.notna(row['Parents']):
            parents = row['Parents'].split(';')
            if len(parents) == 2 and all(p in pos for p in parents):
                p1_x, p1_y = pos[parents[0]]
                p2_x, p2_y = pos[parents[1]]
                mid_x = (p1_x + p2_x) / 2 + 0.5
                ax.plot([p1_x + 0.5, p2_x + 0.5], [p1_y + 1, p2_y + 1], 'k-')
                ax.plot([mid_x, x + 0.5], [p1_y + 1, y + 1], 'k-')

    # Legend
    legend_patches = [Patch(facecolor=color, edgecolor='black', label=label)
                      for label, color in disorder_colors.items()]
    ax.legend(handles=legend_patches, title="Disorders & Status", loc="upper right", fontsize=8, title_fontsize=10)

    # Styling
    ax.set_xlim(-1, max(x for x, y in pos.values()) + 4)
    ax.set_ylim(min(y for x, y in pos.values()) - 3, 2)
    ax.set_aspect('equal')
    ax.axis('off')
    plt.title("Color-Coded Genetic Disorder Pedigree Chart", fontsize=15)
    plt.tight_layout()

    plt.show()

    # --- Disease Summary Report ---
    disorder_summary = {
        'Color Blindness': 0,
        'CB Carrier': 0,
        'Hemophilia': 0,
        'Hemophilia Carrier': 0,
        'Sickle Cell': 0,
        'SC Carrier': 0,
        'Myotonic Dystrophy': 0,
        'Healthy': 0
    }

    for disorders in df['Disorders']:
        for disorder in disorders:
            if disorder in disorder_summary:
                disorder_summary[disorder] += 1

    # --- Print Summary Text Report ---
    print("=== Disease Summary Report ===")
    print(f"Total Individuals: {len(df)}")
    for disorder, count in disorder_summary.items():
        print(f"{disorder}: {count}")
    
    # --- Bar Chart for Summary ---
    labels = list(disorder_summary.keys())
    values = [disorder_summary[k] for k in labels]
    colors = [disorder_colors.get(k, 'gray') for k in labels]

    plt.figure(figsize=(10, 6))
    bars = plt.bar(labels, values, color=colors)
    plt.xticks(rotation=45, ha='right')
    plt.ylabel("Number of Individuals")
    plt.title("Genetic Disorder Summary")
    
    for bar in bars:
        yval = bar.get_height()
        plt.text(bar.get_x() + bar.get_width()/2, yval + 0.3, int(yval), ha='center', fontsize=9)

    plt.tight_layout()
    plt.show()

# Run everything
draw_pedigree(df)
