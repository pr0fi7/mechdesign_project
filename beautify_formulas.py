import json
import re

# Read notebook
with open(r'C:\Users\memor\Desktop\kulproject\machinedesign\lab1\Assignment_1-new.ipynb', 'r', encoding='utf-8') as f:
    nb = json.load(f)

changes = []

# Beautify specific cells
for cell in nb['cells']:
    if cell['cell_type'] != 'markdown':
        continue

    # Get source as string
    if isinstance(cell['source'], list):
        source = ''.join(cell['source'])
    else:
        source = cell['source']

    original = source

    # Beautify patterns - only visual improvements

    # 1. Fix degree symbols
    source = source.replace('0degrees', '0°')
    source = source.replace('degrees', '°')

    # 2. Improve Greek letters in text
    source = source.replace('Beta = 0.0', 'β = 0')
    source = source.replace('beta = 0', 'β = 0')
    source = source.replace('(Beta = 0degrees)', '(β = 0°)')

    # 3. Ensure proper spacing in displayed equations
    # Remove extra spaces before/after $$
    source = re.sub(r'\$\$\s+', '$$\n', source)
    source = re.sub(r'\s+\$\$', '\n$$', source)

    # 4. Standardize inline variable formatting in regular text
    # Only format when NOT already in math mode
    lines = source.split('\n')
    new_lines = []
    for line in lines:
        # Skip if line contains $$ or is already in math
        if '$$' not in line and '$' not in line:
            # Format common variables
            line = re.sub(r'\bK_A\b', r'$K_A$', line)
            line = re.sub(r'\bi_total\b', r'$i_{\text{total}}$', line)
            line = re.sub(r'\bi_pulley\b', r'$i_{\text{pulley}}$', line)
            line = re.sub(r'\bi_gearbox\b', r'$i_{\text{gearbox}}$', line)
            line = re.sub(r'\be\'', r"$e'$", line)
        new_lines.append(line)
    source = '\n'.join(new_lines)

    # 5. Improve multiplication symbol
    source = re.sub(r'(\d+)\s*\*\s*(\d+)', r'\1 × \2', source)

    # 6. Improve range arrows
    source = source.replace(' -> ', ' → ')
    source = source.replace('->', '→')

    # 7. Ensure proper bullet formatting
    source = re.sub(r'^\s*-\s+', '- ', source, flags=re.MULTILINE)

    if source != original:
        # Convert back to list format
        cell['source'] = source.split('\n')
        cell['source'] = [line + '\n' for line in cell['source'][:-1]] + [cell['source'][-1]]
        changes.append(cell['id'][:8])

print(f'Beautified {len(changes)} markdown cells')
if changes:
    print(f'Changed cell IDs: {", ".join(changes)}')

# Save notebook
with open(r'C:\Users\memor\Desktop\kulproject\machinedesign\lab1\Assignment_1-new.ipynb', 'w', encoding='utf-8') as f:
    json.dump(nb, f, indent=1, ensure_ascii=False)

print('Notebook saved successfully')
