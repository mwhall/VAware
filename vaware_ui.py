# -*- coding: utf-8 -*-
import argparse
import random

import dash
import dash_html_components as html
import dash_core_components as dcc
import plotly.graph_objs as go
import plotly
import pandas as pd

def read_table(table_file):
    table_file.seek(0)
    table = pd.read_csv(table_file, sep="\t", header=0, comment="#")
    return table

def read_header(table_file):
    table_file.seek(0)
    line = table_file.readline()
    values = []
    while (line[0] == "#"):
        if (": " in line):
            values.append(line.split(": ")[-1].strip())
        line = table_file.readline()
    return values

def get_unique_taxa(table, level=1):
    tax_set = set()
    tax_list = []
    for full_tax in table["Taxonomy"]:
        tax_split = full_tax.split(";")
        if (len(tax_split) >= level):
            tax = tax_split[level - 1]
            if tax not in tax_set:
                tax_set.add(tax)
                tax_list.append({'label':tax, 'value':tax})
    return tax_list

def generate_random_color():
    r = lambda: random.randint(0,255)
    return '#%02X%02X%02X' % (r(),r(),r())

def update_hist(value):
    fig = plotly.tools.make_subplots(rows=3, cols=3, 
                                     specs=[[{}, {}, {}], [{}, {}, {}], 
                                           [{'colspan': 2}, None, None]],
                                     subplot_titles=['FP Mismatches', 'FP 3\' Mismatches', 
                                                     'FP Gaps', 'RP Mismatches', 
                                                     'RP 3\' Mismatches', 'RP Gaps', 
                                                     'Insert Length'])
    fig['layout']['margin'] = {'l': 30, 'r': 10, 'b': 30, 't': 30}
    fig['layout']['legend'] = {'x': 0, 'y': 1, 'xanchor': 'left'}

#    if col_name not in list(table):
#        raise ValueError("Requested table column {} not found.".format(col_name))

    if type(value) is dict:
        val = [ value['label'] ]
    elif type(value) is list:
        val = value
    else:
        val = [ value ]

    if val is not [ None ]:
        if len(val) > 1:
            fig['layout']['barmode'] = 'overlay'
        fig['layout']['legend'] = {'orientation': 'h'}
        for taxon in val:
            k = 1
            i = 1
            j = 1
            taxon_colour = generate_random_color()
            taxa_filter = [ taxon in x for x in table['Taxonomy'] ]
            subset_table = table.loc[ taxa_filter ]
            for col_name in ['FP Mismatches', 'FP 3\' Mismatches', 'FP Gaps', 'RP Mismatches', 'RP 3\' Mismatches', 'RP Gaps']:
                fig.append_trace({'x': subset_table[col_name].dropna(), 
                                  'type': 'histogram',
                                  'histnorm': 'probability',
                                  'opacity': 0.75,
                                  'marker': {'color': taxon_colour},
                                  'name': taxon,
                                  'showlegend': False},
                                 i, j)
                #k=count, i=rows, j=columns
                k += 1
                i = (k // 4) + 1
                j = (k - 1) % 3 + 1
            fig.append_trace({'x': subset_table['Insert Length'].dropna(),
                              'type': 'histogram',
                              'histnorm': 'probability',
                              'opacity': 0.75,
                              'marker': {'color': taxon_colour},
                              'name': taxon,
                              'showlegend': True},
                             3, 1)

        return fig
    else:
        k=1
        i=1
        j=1
        for col_name in ['FP Mismatches', 'FP 3\' Mismatches', 'FP Gaps', 'RP Mismatches', 'RP 3\' Mismatches', 'RP Gaps']:
            if (i == 2) & (j == 3):
                showlegend = True
            else:
                showlegend = False
            fig.append_trace({'x': table[col_name].dropna(),
                              'type': 'histogram',
                              'opacity': 0.75,
                              'marker': {'color': taxon_colour},
                              'name': taxon,
                              'showlegend': showlegend},
                              i, j)
            k += 1
            i = (k // 4) + 1
            j = (k - 1) % 3 + 1

        return fig


def start_dash(app, input_file):
  params = read_header(input_file)
  tax_list = get_unique_taxa(table)

  app.layout = html.Div([
    html.Div([
      html.Label('Taxonomic Level'),
      dcc.Slider(
        id='taxon-level',
        min=1,
        max=6,
        marks={i: 'Level {}'.format(i) if i == 1 else str(i) for i in range(1, 7)},
        value=1,
      ),
      html.Br(),
      html.Label('Taxonomic Group'),
      dcc.Dropdown(
        id='taxa-list',
        options=tax_list,
        value=tax_list[0],
        multi=True
      ),
      html.Br()
      ], style={'width': '35%'}),
    html.Div([
        dcc.Graph(id='histograms'),
        ])
 ], style={'marginbottom': 50, 'marginTop': 25,
           'marginLeft': 25, 'marginRight': 25})

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("input_file", type=argparse.FileType('r'))
    args = parser.parse_args()
    table = read_table(args.input_file)
    app = dash.Dash()
    start_dash(app, args.input_file)
    
    # **** START REACTIVITY FUNCTIONS ****
    # These go here because they need to be in scope of app
    # but before the server is run
    @app.callback(
      dash.dependencies.Output('taxa-list', 'options'),
      [dash.dependencies.Input('taxon-level', 'value')])
    def update_taxa_list(level):
        return get_unique_taxa(table, level)

    @app.callback(
      dash.dependencies.Output('taxa-list', 'value'),
      [dash.dependencies.Input('taxa-list', 'options')])
    def update_taxa_choice(options):
        return options[0]

    @app.callback(
      dash.dependencies.Output('histograms', 'figure'),
      [dash.dependencies.Input('taxa-list', 'value')])
    def update_histograms(value):
        return update_hist(value)

    app.run_server(debug=True)
