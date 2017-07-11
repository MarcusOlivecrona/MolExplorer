from __future__ import absolute_import

from bokeh.core.properties import Bool, Color, Either, Enum, Float, Instance, Int, List, Override, String
from bokeh.models.widgets.tables import TableWidget, TableColumn

class HighTable(TableWidget):
    ''' Two dimensional grid for visualisation and editing large amounts
    of data.
    '''

    __implementation__ = "high_table.coffee"

    columns = List(Instance(TableColumn), help="""
    The list of child column widgets.
    """)

    fit_columns = Bool(True, help="""
    Whether columns should be fit to the available width. This results in no
    horizontal scrollbar showing up, but data can get unreadable if there is
    no enough space available. If set to ``True``, columns' width is
    understood as maximum width.
    """)

    sortable = Bool(True, help="""
    Allows to sort table's contents. By default natural order is preserved.
    To sort a column, click on it's header. Clicking one more time changes
    sort direction. Use Ctrl + click to return to natural order. Use
    Shift + click to sort multiple columns simultaneously.
    """)

    reorderable = Bool(True, help="""
    Allows the reordering of a tables's columns. To reorder a column,
    click and drag a table's header to the desired location in the table.
    The columns on either side will remain in their previous order.
    """)

    editable = Bool(False, help="""
    Allows to edit table's contents. Needs cell editors to be configured on
    columns that are required to be editable.
    """)

    selectable = Either(Bool(True), Enum("checkbox"), help="""
    Whether a table's rows can be selected or not. Using ``checkbox`` is
    equivalent  to ``True``, but makes selection visible through a checkbox
    for each row,  instead of highlighting rows. Multiple selection is
    allowed and can be achieved by either clicking multiple checkboxes (if
    enabled) or using Shift + click on rows.
    """)

    row_headers = Bool(True, help="""
    Enable or disable row headers, i.e. the index column.
    """)

    scroll_to_selection = Bool(True, help="""
    Whenever a selection is made on the data source, scroll the selected
    rows into the table's viewport if none of the selected rows are already
    in the viewport.
    """)

    height = Override(default=400)
