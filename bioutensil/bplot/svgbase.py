import svgwrite


class SVGObject():
    def __init__(self, insert, size,
                 anchor_object=None, anchor_pos='',
                 target_object=None, target_pos='',
                 **kwargs
                 ):

        self.insert = insert
        self.size = size

        self.anchor_object = anchor_object
        self.anchor_pos = anchor_pos
        self.target_object = target_object
        self.target_pos = target_pos

        if self.anchor_object:
            self._update_insert()
        if self.target_object:
            self._update_size()

    def _update_insert(self):
        x, y = self.anchor_object.x, self.anchor_object.y
        width, height = self.anchor_object.width, self.anchor_object.height

        if 'iup' in self.anchor_pos:
            self.y = y
        if 'oup' in self.anchor_pos:
            self.y = y - self.height

        if 'center' in self.anchor_pos:
            self.y = y + height/2 - self.height/2

        if 'ibottom' in self.anchor_pos:
            self.y = y + height - self.height
        if 'obottom' in self.anchor_pos:
            self.y = y + height

        if 'ileft' in self.anchor_pos:
            self.x = x
        if 'oleft' in self.anchor_pos:
            self.x = x - self.width

        if 'middle' in self.anchor_pos:
            self.x = x + width/2 - self.width/2

        if 'iright' in self.anchor_pos:
            self.x = x + width - self.width
        if 'oright' in self.anchor_pos:
            self.x = x + width

        self.insert = (self.x, self.y)

    def _update_size(self):
        x, y = self.target_object.x, self.target_object.y
        width, height = self.target_object.width, self.target_object.height

        if 'iup' in self.target_pos:
            self.y2 = y
        if 'oup' in self.target_pos:
            self.y2 = y - self.height

        if 'center' in self.target_pos:
            self.y2 = y + height/2 - self.height/2

        if 'ibottom' in self.target_pos:
            self.y2 = y + height - self.height
        if 'obottom' in self.target_pos:
            self.y2 = y + height

        if 'ileft' in self.target_pos:
            self.x2 = x
        if 'oleft' in self.target_pos:
            self.x2 = x - self.width

        if 'middle' in self.target_pos:
            self.x2 = x + width/2 - self.width/2

        if 'iright' in self.target_pos:
            self.x2 = x + width - self.width
        if 'oright' in self.target_pos:
            self.x2 = x + width

        self.size = (self.x2 - self.x, self.y2 - self.y)


class Rect(SVGObject):
    def __init__(self, insert, size,
                 anchor_object=None, anchor_pos='',
                 target_object=None, target_pos='',
                 **kwargs):

        self.x, self.y = insert
        self.width, self.height = size

        super(Rect, self).__init__(insert, size,
                                   anchor_object=anchor_object,
                                   anchor_pos=anchor_pos,
                                   target_object=target_object,
                                   target_pos=target_pos,
                                   **kwargs)
        self.svgobj = svgwrite.shapes.Rect(self.insert, self.size, **kwargs)


class Circle(SVGObject):
    def __init__(self, insert, size,
                 anchor_object=None, anchor_pos='',
                 target_object=None, target_pos='',
                 **kwargs):

        self.x, self.y = insert
        self.width, self.height = size

        super(Circle, self).__init__(insert, size,
                                     anchor_object=anchor_object,
                                     anchor_pos=anchor_pos,
                                     target_object=target_object,
                                     target_pos=target_pos,
                                     **kwargs)

        self.r = size[0]/2.0
        center = self.insert[0] + self.r, self.insert[1] + self.r
        self.svgobj = svgwrite.shapes.Circle(center, self.r, **kwargs)


class Line(SVGObject):
    def __init__(self, insert, size,
                 anchor_object=None, anchor_pos='',
                 target_object=None, target_pos='',
                 **kwargs):
        self.x, self.y = insert
        self.width, self.height = size

        super(Line, self).__init__(insert, size,
                                   anchor_object=anchor_object,
                                   anchor_pos=anchor_pos,
                                   target_object=target_object,
                                   target_pos=target_pos,
                                   **kwargs)

        width, height = self.size
        end = self.x + width, self.y + height
        self.svgobj = svgwrite.shapes.Line(self.insert, end, **kwargs)


class Text(SVGObject):

    def __init__(self, text, insert=(0, 0),
                 font_size=10, font_family='Arial', font_style='normal',
                 anchor_object=None, anchor_pos='',
                 target_object=None, target_pos='',
                 **kwargs):

        self.x, self.y = insert
        self.width, self.height = 0, 0

        super(Text, self).__init__(insert, None,
                                   font_size=font_size,
                                   font_family=font_family,
                                   font_style=font_style,
                                   anchor_object=anchor_object,
                                   anchor_pos=anchor_pos,
                                   target_object=target_object,
                                   target_pos=target_pos,
                                   **kwargs)

        self.svgobj = svgwrite.text.Text(text, self.insert,
                                         font_size=font_size,
                                         font_family=font_family,
                                         font_style=font_style,
                                         dy=self.dy, style=self.style,
                                         **kwargs)

    def _update_insert(self):
        # font render in up-right direction
        SVGObject._update_insert(self)
        self.dy = ['1em']
        self.style = 'text-anchor: start'

        if 'iup' in self.anchor_pos:
            self.dy = ['1em']
        if 'oup' in self.anchor_pos:
            self.dy = ['-0.3em']

        if 'center' in self.anchor_pos:
            self.dy = ['0.3em']

        if 'ibottom' in self.anchor_pos:
            self.dy = ['-0.3em']
        if 'obottom' in self.anchor_pos:
            self.dy = ['1em']

        offset = 2
        if 'ileft' in self.anchor_pos:
            self.style = 'text-anchor: start'
            self.x += offset
        if 'oleft' in self.anchor_pos:
            self.style = 'text-anchor: end'
            self.x -= offset

        if 'iright' in self.anchor_pos:
            self.style = 'text-anchor: end'
            self.x -= offset
        if 'oright' in self.anchor_pos:
            self.style = 'text-anchor: start'
            self.x += offset

        if 'middle' in self.anchor_pos:
            self.style = 'text-anchor: middle'

        self.insert = self.x, self.y


class LinearGradientRect(Rect):

    def __init__(self, insert, size, start_color, start_opacity, end_color, end_opacity, **extra):
        super(LinearGradientRect, self).__init__(insert, size, **extra)

        self._start_color = start_color
        self._start_opacity = start_opacity
        self._end_color = end_color
        self._end_opacity = end_opacity

    def render(self):
        horizontal_gradient = self._dwg.linearGradient((0, 0), (1, 0))
        horizontal_gradient.add_stop_color(0, self._start_color, self._start_opacity)
        horizontal_gradient.add_stop_color(1, self._end_color, self._end_opacity)

        rect = self._dwg.rect(self.insert, self.size, fill=horizontal_gradient.get_paint_server(default='currentColor'))

        return horizontal_gradient, rect


class Container(SVGObject):
    def __init__(self, insert, size,
                 h_margin=0, v_margin=0,
                 caption=None,
                 caption_pos='',
                 fill='none',
                 border_stroke='none',
                 stroke='none',
                 anchor_object=None,
                 anchor_pos='',
                 target_object=None,
                 target_pos='',
                 **kwargs):

        self.h_margin = h_margin
        self.v_margin = v_margin
        width, height = size
        size = width+2*h_margin, height+2*v_margin

        self.caption = caption
        self.caption_pos = caption_pos
        self.x, self.y = size
        self.width, self.height = size

        super(Container, self).__init__(insert, size,
                                        anchor_object=anchor_object,
                                        anchor_pos=anchor_pos,
                                        target_object=target_object,
                                        target_pos=target_pos,
                                        **kwargs)
        self.x, self.y = self.size
        self.width, self.height = self.size

        self.svgobj = svgwrite.container.SVG(self.insert, self.size,
                                             stroke='none',
                                             **kwargs)  # update the data
        self.rect = Rect((0, 0), self.size,
                         fill=fill,
                         stroke=border_stroke,
                         **kwargs)  # update the data
        self.svgobj.add(self.rect.svgobj)

        if caption:
            self._draw_caption()

    def _draw_caption(self):
        label = Text(self.caption, (0, 0),
                     font_size=8, font_weight='bold',
                     anchor_object=self.rect,
                     anchor_pos=self.caption_pos)
        self.svgobj.add(label.svgobj)

    def split(self, split_pattern=None, split_ratio=None):
        width, height = self.size
        if split_pattern:
            split_x, split_y = split_pattern
            split_ratio_x = [1.0/split_x] * split_x
            split_ratio_y = [1.0/split_y] * split_y
        if split_ratio:
            split_ratio_x, split_ratio_y = split_ratio
            split_ratio_x = list(map(int, split_ratio_x.split(':')))
            split_ratio_y = list(map(int, split_ratio_y.split(':')))
            split_ratio_x = list(map(lambda n: 1.0*n/sum(split_ratio_x), split_ratio_x))
            split_ratio_y = list(map(lambda n: 1.0*n/sum(split_ratio_y), split_ratio_y))

        items = []
        child_y = 0
        for j in range(len(split_ratio_y)):
            child_height = height*split_ratio_y[j]

            line = []
            child_x = 0
            for i in range(len(split_ratio_x)):
                child_width = width*split_ratio_x[i]

                insert = child_x, child_y
                size = child_width, child_height
                child = Container(insert, size)
                line.append(child)

                child_x += child_width

            items.append(line)
            child_y += child_height
        return items
