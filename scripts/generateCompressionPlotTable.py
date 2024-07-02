import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns

if (True):
    # Data for the table
    columns = ("Species", "Number of\nSequences", "Average\nGenome Length", "Average Mash\nDistance" ,"PanMAN\nFile Size (MB)", "GFA", "VG",  "GBZ","Pangraph")
    data = [
        ["Sars-CoV-2", "20,000", "29,902", "0.0025","2.0","523.3X", "153.5X",  "17.1X", "614.4X"],
        ["HIV", "20,000", "9,299", "0.091","9.1", "135.1X", "42.2X",  "6.8X","49.6X"],
        [ "RSV", "4,000", "14,987", "0.123", "0.5","680.4X", "221.7X",  "5.2X", "234.8X"],
        ["E. Coli", "1,000", "4,641,650", "0.023","114.0","21.6X", "6.6X",  "7.9X", "30.5X"],
        [ "M. Tuberculosis", "400", "4,411,448",  "0.00093","5.0", "148.8X", "79.6X", "42.1X","80.4X"],
        [ "K. pneumoniae", "1,000", "5,328,875",  "0.023","200.0", "19.5X", "6.1X", "7.7X", "26.1X"]
        
    ]

    # Create figure for the table
    fig, ax = plt.subplots(figsize=(16, 3))
    ax.axis('tight')
    ax.axis('off')

    # Create the table and adjust layout
    table = ax.table(cellText=data, colLabels=columns, cellLoc = 'center', loc='center', colWidths=[0.14,0.09,0.12,0.1,0.1, 0.08, 0.08, 0.08, 0.08, 0.09])
    table.auto_set_font_size(False)
    table.set_fontsize(12)
    table.scale(1.2, 1.4)

    # Adjust the height of the first row cells
    first_row_height = 0.3  # Increase this value to increase the height of the first row
    for pos, cell in table.get_celld().items():
        if pos[0] == 0:  # Check if it is the first row (header row)
            cell.set_height(first_row_height)
        # elif pos[0] == 6:
        #     cell.set_height(0.2)
        elif pos[0] in [1,2,3,4,5,6]:
            cell.set_height(0.1)


    # Bold the first column and first row
    for cell in table.get_celld().values():
        cell.set_text_props(fontweight='normal')  # Resetting all to normal before setting first column to bold
    for (i, j), cell in table.get_celld().items():
        if i == 0 or i == -1:
            cell.set_text_props(fontweight='bold')
        if j == 0 or j == -1:
            cell.set_text_props(fontweight='bold')

    # Adding title
    # ax.set_title("Compression achieved by PanMAN compared to other formats ", fontsize=16, weight='bold', pad=5)

    plt.savefig("../figures/compressionTable.svg", format="svg")
