import svgwrite


class SVGObject(object):
    def __init__(self, insert, size,
                 fill='white', opacity=1,
                 stroke='black', stroke_width=1, stroke_dasharray=[], stroke_lines=['outline'],
                 parent=None,
                 anchor_object=None,
                 related_insert='',
                 self_pos='',
                 rotate=None,
                 transform=None,
                 profile='full',
                 **kwargs
                 ):

        self._dwg = svgwrite.Drawing(profile=profile)

        self.insert = insert
        self.x, self.y = insert
        self.size = size

        self.fill = fill
        self.opacity = opacity

        self.stroke = stroke
        self.stroke_width = stroke_width
        self._stroke_lines = stroke_lines
        if stroke_dasharray:
            self.stroke_dasharray = stroke_dasharray
        elif stroke_lines and self.size:
            self._draw_stroke_lines()

        # text
        self.dx = None
        self.dy = None
        self.style = None

        # line

        self.rotate = rotate
        self.transform = transform
        if rotate is not None:
            self._update_rotate()

        # the last part
        self.parent = parent
        self.anchor_object = anchor_object
        self.related_insert = related_insert
        self.self_pos = self_pos

        self._update_insert()
        self._update_pos()

    def _update_insert(self):
        # use anchor related position

        if self.anchor_object:  # absolute position
            if isinstance(self.anchor_object, svgwrite.shapes.Rect):
                x, y = self.anchor_object['x'], self.anchor_object['y']
                width, height = self.anchor_object['width'], self.anchor_object['height']
            if isinstance(self.anchor_object, svgwrite.shapes.Line):
                # line object
                x, y = self.anchor_object['x1'], self.anchor_object['y1']
                x2, y2 = self.anchor_object['x2'], self.anchor_object['y2']
                width, height = x2 - x, y2 - y
            if isinstance(self.anchor_object, svgwrite.shapes.Circle):
                r = self.anchor_object['r']
                height, width = r * 2, r * 2
                x, y = self.anchor_object['cx'] - r, self.anchor_object['cy'] - r

        if self.parent: # related position based on parent
            x, y = 0, 0
            width, height = self.parent['width'], self.parent['height']

        if 'up' in self.related_insert:             self.y = y
        elif 'center' in self.related_insert:       self.y = y + height/2.0
        elif 'bottom' in self.related_insert:       self.y = y + height
        elif self.anchor_object:                self.y += y

        if 'left' in self.related_insert:           self.x = x
        elif 'middle' in self.related_insert:       self.x = x + width/2.0
        elif 'right' in self.related_insert:        self.x = x + width
        elif self.anchor_object:                self.x += x

        self.insert = (self.x, self.y)

    def _update_pos(self):
        if self.self_pos:
            self._update_pos_by_self_pos()
        else:
            pass
            # self._update_pos_by_related_insert()

    def _update_pos_by_related_insert(self):
        width, height = self.size
        if 'middle' in self.related_insert:
            self.x -= width/2.0
        if 'right' in self.related_insert:
            self.x -= width
        if 'center' in self.related_insert:
            self.y -= height/2.0
        if 'bottom' in self.related_insert:
            self.y -= height

        self.insert = (self.x, self.y)

    def _update_pos_by_self_pos(self):
        width, height = self.size
        if 'right' in self.self_pos:      self.x
        if 'middle' in self.self_pos:     self.x -= width/2.0
        if 'left' in self.self_pos:      self.x -= width
        if 'bottom' in self.self_pos:     self.y
        if 'center' in self.self_pos:     self.y -= height/2.0
        if 'up' in self.self_pos:     self.y -= height

        self.insert = (self.x, self.y)

    def _update_rotate(self):
        # horizontal as default
        if self.rotate is not None:
            transform = 'rotate({} {} {})'.format(self.rotate, self.x, self.y)
            self.transform = transform

    def _draw_stroke_lines(self):
        width, height = self.size
        if 'outline' in self._stroke_lines: self._stroke_lines += ['up','right','bottom','left']
        pattern = [int(x in self._stroke_lines) for x in ['up','right','bottom','left']]
        size = [width, height, width, height]
        stroke_dasharray = []
        prev_pattern = -1
        element = 0
        for pattern, size in zip(pattern, size):
            if prev_pattern == pattern or prev_pattern == -1:
                element += size
            else: # new
                # update prev
                stroke_dasharray.append(element)
                element = size
            prev_pattern = pattern


        stroke_dasharray.append(element)

        if 'up' not in self._stroke_lines:
            stroke_dasharray = [0] + stroke_dasharray

        self.stroke_dasharray=','.join(map(str,stroke_dasharray))


class Rect(SVGObject):
    def __init__(self, insert, size,
            rx=0, ry=0, fill='white', opacity=1,
            **extra):
        super(Rect, self).__init__(insert, size, fill=fill, opacity=opacity, **extra)
        self.rx = rx
        self.ry = ry


    def render(self):
        rect = self._dwg.rect(self.insert, self.size,
                rx=self.rx, ry=self.ry,
                stroke=self.stroke, stroke_width=self.stroke_width, stroke_dasharray=self.stroke_dasharray,
                fill=self.fill, opacity=self.opacity)
        return rect


class Circle(SVGObject):
    def __init__(self, insert, size,
                 fill='white', opacity=1,
                 **extra):
        super(Circle, self).__init__(insert, size, fill=fill, opacity=opacity, **extra)
        self.r = size[0]/2.0

    def render(self):
        center = self.insert[0] + self.r, self.insert[1] + self.r
        rect = self._dwg.circle(center, self.r,
                                stroke=self.stroke, stroke_width=self.stroke_width, stroke_dasharray=self.stroke_dasharray,
                                fill=self.fill, opacity=self.opacity)
        return rect


class Line(SVGObject):
    def __init__(self, insert, size, **extra):
        super(Line, self).__init__(insert, size, **extra)

    def render(self):
        width, height = self.size
        end = self.x + width, self.y + height
        return self._dwg.line(self.insert, end, stroke=self.stroke)


class Text(SVGObject):

    def __init__(self, text, insert=(0, 0),
                 font_size=10, font_family='Arial', font_style='normal',
                 fill='black', letter_spacing=0,
                 **extra):
        super(Text, self).__init__(insert, None, **extra)
        self.text = text
        self.font_size = font_size
        self.font_family = font_family
        self.font_style = font_style
        self.fill = fill
        self.letter_spacing = letter_spacing

    def _update_pos(self):
        if 'up' in self.related_insert:
            self.dy = ['-0.3em']
        if 'center' in self.related_insert:
            self.dy = ['0.3em']
        if 'bottom' in self.related_insert:
            self.dy = ['1em']

        if 'left' in self.related_insert:
            self.style = 'text-anchor: end'
        if 'right' in self.related_insert:
            self.style = 'text-anchor: start'
        if 'middle' in self.related_insert:
            self.style = 'text-anchor: middle'
        if 'center' not in self.related_insert:
            self.style = 'text-anchor: middle'

        if self.style == 'text-anchor: end':
            self.dx = ['-0.5em']
        if self.style == 'text-anchor: start':
            self.dx = ['0.5em']

    def render(self):
        return self._dwg.text(self.text, self.insert, dx=self.dx, dy=self.dy,
                              style=self.style,
                              fill=self.fill,
                              font_size=self.font_size,
                              font_family=self.font_family,
                              font_style=self.font_style,
                              letter_spacing=self.letter_spacing)




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
    def __init__(self, insert=(0,0), size=(10,10), stroke='gray',
            h_margin=0, v_margin=0,
            **extra):
        self.stroke = stroke

        self.h_margin = h_margin
        self.v_margin = v_margin
        width, height = size
        size = width+2*h_margin, height+2*v_margin

        super(Container, self).__init__(insert, size, **extra)
        self._svg = svgwrite.container.SVG(self.insert, self.size) # update the data


        #self.draw_lines('right','left','up','bottom')

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

    def draw_lines(self, *names):
        names = list(names)
        if 'outline' in names:
            names += ['left', 'right', 'up', 'bottom']
        width, height = self.size
        if 'left' in names:
            line = Line((0,0), (0,height), parent=self._svg, stroke=self.stroke).render()
            self._svg.add(line)
        if 'right' in names:
            line = Line((width,0), (0,height), parent=self._svg, stroke=self.stroke).render()
            self._svg.add(line)
        if 'middle' in names:
            line = Line((width/2,0), (0,height), parent=self._svg, stroke=self.stroke).render()
            self._svg.add(line)

        if 'up' in names:
            line = Line((0,0), (width,0), parent=self._svg, stroke=self.stroke).render()
            self._svg.add(line)
        if 'bottom' in names:
            line = Line((0,height), (width,0), parent=self._svg, stroke=self.stroke).render()
            self._svg.add(line)
        if 'center' in names:
            line = Line((0,height/2), (width,0), parent=self._svg, stroke=self.stroke).render()
            self._svg.add(line)

    def render(self):
        return self._svg




