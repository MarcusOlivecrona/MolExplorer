from bokeh.plotting import figure, ColumnDataSource, curdoc
from bokeh.models import HoverTool, PanTool, BoxZoomTool, ResetTool, BoxSelectTool, LassoSelectTool, TapTool, CustomJS
from bokeh.models.widgets import DataTable, TableColumn, HTMLTemplateFormatter, StringEditor, Button
from high_table import HighTable
from bokeh.layouts import row, column, widgetbox, layout
from bokeh.models.widgets import Select
import bokeh.palettes
from rdkit import Chem
from rdkit.Chem import Draw
import sys
import os

error_msg = """Need to provide valid SMILES file as first argument.
                     'bokeh serve . --args [smiles_file]'"""
try:
    path = sys.argv[1]
except IndexError:
    raise IndexError(error_msg)
if not os.path.isfile(path):
    raise ValueError(error_msg)

with open(path, "r") as f:
    smiles = []
    imgs = []

    # For the first header line
    headers = f.readline().split()
    properties = {header:list() for header in headers[1:]}

    print("Property headers: {}".format([h for h in properties.keys()]))

    for line in f:
        try:
            line = line.split()
            mol = Chem.MolFromSmiles(line[0])
            if mol:
                if len(line[0])>70:
                    smiles.append("{:.67}...".format(line[0]))
                else:
                    smiles.append(line[0])
                imgs.append(Draw.MolsToGridImage((mol,), molsPerRow=1, subImgSize=(300,300), useSVG=True))
                for i, property in enumerate(properties.values()):
                    property.append(float(line[i + 1]))
        except IndexError:
            continue

    if len(properties)==0:
        import numpy as np
        print("File did not have any property columns, adding dummy properties")
        n = int(np.sqrt(len(smiles)))
        properties[' '] = [i % n for i, _ in enumerate(smiles)]
        properties['  '] = [i // n for i, _ in enumerate(smiles)]
    elif len(properties)==1:
        print("File only had one property column, adding a dummy property")
        properties[' '] = [i for i, _ in enumerate(smiles)]

    table_data = {key : item for key, item in properties.items()}
    table_data['SMILES'] = smiles
    table_data['IMGS'] = imgs
    print("Reading data completed...")

hover = HoverTool(tooltips="""
    <div>
        <div>
            @IMGS{safe}
        </div>
        <div>
            <span style="font-size: 12px; font-style: bold; color: #404040;">@SMILES</span>
        </div>
        <div>
            <span style="font-size: 15px; color: #F17022;">@x_desc : @x</span>
        </div>
        <div>
            <span style="font-size: 15px; color: #F17022;">@y_desc : @y</span>
        </div>
        <div>
            <span style="font-size: 15px; color: #F17022;">@color_desc : @color_data</span>
        </div>
    </div>
    """, attachment="horizontal")

source = ColumnDataSource(data=dict(x=[], y=[], color=[], color_data=[], x_desc=[], y_desc=[],
                                    color_desc=[], SMILES=[], IMGS=[]))

x_axis = Select(title="X Axis", options=list(properties.keys()), value=list(properties.keys())[0])
y_axis = Select(title="Y Axis", options=list(properties.keys()), value=list(properties.keys())[1])
color_axis = Select(title="Color using", options=list(properties.keys()) + ["None"], value="None")

p = figure(plot_width=1000, plot_height=800, tools=[hover, PanTool(), BoxZoomTool(), TapTool(),
                                                    BoxSelectTool(), LassoSelectTool(), ResetTool()])
p.circle("x", "y", size=20, line_color="color", fill_color="color", fill_alpha=0.5, source=source)

formatter =  HTMLTemplateFormatter(template="<div><%= value %></div>")
table_source = ColumnDataSource(data=table_data)
columns = [TableColumn(field='IMGS', title='Structure', width=300, formatter=formatter),
           TableColumn(field='SMILES', title='SMILES', width=300, editor=StringEditor())] + \
          [TableColumn(field=prop, title=prop, width=150) for prop in properties.keys()]
table = HighTable(source=table_source, columns=columns, editable=True, width=1200, height=1200)

button = Button(label="Download", button_type="warning")
button.callback = CustomJS(args=dict(source=table_source),
                           code=open(os.path.join(os.path.dirname(__file__), "download.js")).read())

def colors_from_data(data):
    min_data = min(data)
    max_data = max(data)
    if not max_data == min_data:
        data_0_to_256 = [int(255 * (data_point - min_data) / (max_data - min_data)) for data_point in data]
        colors = [bokeh.palettes.Plasma256[i] for i in data_0_to_256]
    else:
        colors = [bokeh.palettes.Plasma256[200] for _ in data]
    return colors

def update():
    x_data = properties[x_axis.value]
    y_data = properties[y_axis.value]

    if color_axis.value == "None":
        color_data = [0] * len(x_data)
    else:
        color_data = properties[color_axis.value]

    color = colors_from_data(color_data)

    p.xaxis.axis_label = x_axis.value
    p.yaxis.axis_label = y_axis.value

    x_desc = [x_axis.value] * len(x_data)
    y_desc = [y_axis.value] * len(x_data)
    color_desc = [color_axis.value] * len(x_data)

    source.data = dict(x=x_data, y=y_data, color=color, color_data=color_data, x_desc=x_desc,
                       y_desc=y_desc, color_desc=color_desc, SMILES=smiles, IMGS=imgs)

def update_table(attr, old, new):
    idxs = new['1d']['indices']
    data = {key : [item[i] for i in idxs] for key, item in table_data.items()}
    table_source.data = data

x_axis.on_change("value", lambda attr, old, new: update())
y_axis.on_change("value", lambda attr, old, new: update())
color_axis.on_change("value", lambda attr, old, new: update())
source.on_change("selected", update_table)

widgetbox = widgetbox([x_axis, y_axis, color_axis], width=200)
layout = layout([[p, widgetbox], [table], [button]])

update()
curdoc().add_root(layout)
