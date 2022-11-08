import svgwrite
from bplot import svgbase
import math


class SVGObject(object):

    def __init__(self,
            stroke = 'black', stroke_width = 1,
            stroke_dasharray=[], stroke_lines=['outline'], stroke_pattern=[1,1,1,1],
            fill='white', opacity=1,
            ):
        self._dwg = svgwrite.Drawing()

        self.stroke = stroke
        self.stroke_width = stroke_width
        self._stroke_lines = stroke_lines
        self._stroke_pattern = stroke_pattern
        if stroke_lines or stroke_pattern:
            self._draw_stroke_lines()

        self.fill = fill
        self.opacity = opacity


class Polygon(SVGObject):
    def __init__(self, points, **extra):
        self.points = points
        super(Polygon, self).__init__(**extra)


    def render(self):
        polygon = self._dwg.polygon(self.points,
                stroke=self.stroke, stroke_width=self.stroke_width, stroke_dasharray=self.stroke_dasharray,
                fill=self.fill, opacity=self.opacity)
        return polygon

    def _draw_stroke_lines(self):
        if 'outline' in self._stroke_lines: self._stroke_lines += ['up','right','bottom','left']
        pattern = [int(x in self._stroke_lines) for x in ['up','right','bottom','left']]
        if self._stroke_pattern: pattern = self._stroke_pattern

        distance = lambda x, y: math.sqrt(math.pow(x[0]-y[0], 2) + math.pow(x[1]-y[1], 2))
        size = [distance(x,y) for x, y in zip(self.points, self.points[1:]+[self.points[0]])]
        stroke_dasharray = []
        prev_pattern = -1
        element = 0
        for x, size in zip(pattern, size):
            if prev_pattern == x or prev_pattern == -1:
                element += size
            else: # new
                # update prev
                stroke_dasharray.append(element)
                element = size
            prev_pattern = x

        stroke_dasharray.append(element)

        if pattern[0] != 1:
            stroke_dasharray = [0] + stroke_dasharray

        self.stroke_dasharray=','.join(map(str,stroke_dasharray))


class Arrow(svgbase.Container):

    def __init__(self, insert, size,
                 border_stroke='none',
                 anchor_object=None,
                 anchor_pos='',
                 **kwargs):
        super(Arrow, self).__init__(insert, size,
                                    anchor_object=anchor_object,
                                    anchor_pos=anchor_pos,
                                    border_stroke=border_stroke,
                                    **kwargs)

        line = svgbase.Line((0, 0), (0, 0),
                            anchor_object=self.rect,
                            anchor_pos='middle:iup',
                            target_object=self.rect,
                            target_pos='middle:ibottom',
                            **kwargs)
        self.svgobj.add(line.svgobj)

        width = self.width/2
        line = svgbase.Line((0, 0), (width, width),
                            anchor_object=self.rect,
                            anchor_pos='ileft:ibottom',
                            **kwargs)
        self.svgobj.add(line.svgobj)

        line = svgbase.Line((0, 0), (width, width),
                            anchor_object=self.rect,
                            anchor_pos='iright:ibottom',
                            **kwargs)
        self.svgobj.add(line.svgobj)

