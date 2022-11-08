import svgwrite

import svgobject

class ColorBar(svgobject.Container):

    def __init__(self, insert, size, name, min_value, max_value, color, font_size=10, **extra):

        self.name = name
        self.min_value = min_value
        self.max_value = max_value
        self.color = color
        self.font_size = font_size

        super(ColorBar, self).__init__(insert, size, **extra)

    def render(self):

        size = (50,20)
        related_insert = ('center', 'middle')
        defs, colorbar = svgobject.LinearGradientRect(self.insert, size, self.color, 0.3, self.color, 1.0, parent=self._svg, related_insert=related_insert).render()
        self._svg.add(colorbar)
        self._svg.defs.add(defs)

        related_insert = ('up', 'middle')
        name_text = svgobject.Text(self.name, self.insert, font_size=self.font_size, anchor_object=colorbar, related_insert=related_insert).render()
        self._svg.add(name_text)

        min_value_str = str(self.min_value)
        related_insert = ('bottom', 'left')
        min_value_text = svgobject.Text(min_value_str, self.insert, font_size=self.font_size, anchor_object=colorbar, related_insert=related_insert).render()
        self._svg.add(min_value_text)


        max_value_str = str(round(self.max_value, 4))
        related_insert = ('bottom', 'right')
        max_value_text = svgobject.Text(max_value_str, self.insert, font_size=self.font_size, anchor_object=colorbar, related_insert=related_insert).render()
        self._svg.add(max_value_text)

        return self._svg





