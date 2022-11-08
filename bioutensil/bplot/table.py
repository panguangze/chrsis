from bplot import svgbase
from bplot import text


class Table(svgbase.Container):

    def __init__(self, insert, size, matrix,
                 rows=None,
                 cols=None,
                 show_value=None,
                 show_rows=None,
                 show_cols=None,
                 value_pos='middle:center',
                 show_grid=None,
                 font_size=10,
                 cell_fill='none',
                 value_fill='black',
                 value_cmap=None,
                 cell_cmap=None,
                 **kwargs
                 ):
        self.matrix = matrix
        self.rows = rows
        self.cols = cols
        self.show_value = show_value
        self.value_pos = value_pos
        self.show_rows = show_rows
        self.show_cols = show_cols

        self.show_grid = show_grid

        self.font_size = font_size

        self.value_fill = value_fill
        self.cell_fill = cell_fill
        self.value_cmap = value_cmap or []
        self.cell_cmap = cell_cmap or []

        self.ratio_x = ':'.join(['1']*len(self.matrix[0]))
        self.ratio_y = ':'.join(['1']*len(self.matrix))
        super(Table, self).__init__(insert, size,
                                    **kwargs)

        self._render()
        if rows:
            self._render_rows()
        if cols:
            self._render_cols()

    def _render(self):

        rows = self.split(split_ratio=(self.ratio_x, self.ratio_y))
        for i, row in enumerate(rows):
            for j, cell in enumerate(row):
                if self.show_grid:
                    cell_stroke = 'black'
                else:
                    cell_stroke = 'white'

                if self.cell_cmap:
                    cell_fill = self.cell_cmap[i][j]
                else:
                    cell_fill = self.cell_fill

                rect = svgbase.Rect((0, 0), cell.size,
                                    fill=cell_fill, stroke=cell_stroke)
                cell.svgobj.add(rect.svgobj)
                cell.rect = rect

                v = self.matrix[i][j]
                if self.value_cmap:
                    value_fill = self.value_cmap[i][j]
                else:
                    value_fill = self.value_fill
                if self.show_value:
                    self._set_label(cell, v,
                                    show_value=self.show_value,
                                    value_pos=self.value_pos,
                                    fill=value_fill)
                self.svgobj.add(cell.svgobj)

    def _render_rows(self):
        width = 10
        for row in self.rows:
            w = text.measure_text(''.join(str(x) for x in row), self.font_size) + 2*len(row)
            if w > width:
                width = w

        x, y = self.insert
        _, height = self.size
        insert = x - width, y
        size = width, height

        container = Table(insert, size, self.rows,
                          show_value=self.show_rows,
                          value_pos='iright:center',
                          font_size=self.font_size,
                          )
        self.trows = container

    def _render_cols(self):
        height = self.font_size*2

        x, y = self.insert
        width, _ = self.size
        insert = x, y - height
        size = width, height
        container = Table(insert, size, self.cols,
                          show_value=self.show_cols,
                          value_pos='middle:ibottom',
                          font_size=self.font_size,
                          )
        self.tcols = container

    def _set_label(self, cell, v,
                   show_value=None,
                   value_pos='',
                   fill='black',
                   ):
        try:
            v = round(v, 2)
        except Exception:
            pass
        if show_value == 'latex':
            name = v.replace('\\', '').replace('{', '').replace('}', '')
            size = cell.size[0]*2/3, cell.size[1]*2/3
            label = text.LatexText(v, name, (0, 0), size,
                                   font_size=self.font_size, fill='black',
                                   anchor_object=cell.rect, anchor_pos=value_pos,
                                   )
            cell.svgobj.add(label.svgobj)
            return

        if show_value == 'sign':
            if v > 0:
                v = '+'
            elif v < 0:
                v = '-'
            else:
                v = '0'

        label = svgbase.Text(v,
                             anchor_object=cell.rect,
                             anchor_pos=value_pos,
                             fill=fill,
                             font_size=self.font_size,
                             )
        cell.svgobj.add(label.svgobj)
