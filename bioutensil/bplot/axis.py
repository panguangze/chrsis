import math


import svgobject


class Axis(svgobject.Container):

    def __init__(self, insert, width, ticks_gap_width, 
            label='axis label', label_width=15, tick_self_width=3, height=10,
            mirror=False, h_margin=80, v_margin=10,
            **extra):
        self.ticks_gap_width = ticks_gap_width
        self.tick_self_width = tick_self_width

        self.tick_num = int(math.ceil(width/ticks_gap_width)) + 1
        self.length =  (self.tick_num - 1) * ticks_gap_width
        self.label = label
        self.label_width = label_width
        self.mirror = mirror

        size = self.length, height
        super(Axis, self).__init__(insert, size, 
                h_margin=h_margin, v_margin=v_margin, **extra)


    def render(self):
        related_insert = ('center', 'middle')
        insert = (0, 0)
        size = (self.length, 0)
        line = svgobject.Line(insert, size, stroke='black', parent=self._svg, related_insert=related_insert).render()
        # add axis label
        related_insert = ('center', 'right')
        label_text = svgobject.Text(self.label, self.insert, anchor_object=line, related_insert=related_insert).render()

        for i in range(self.tick_num):
            _, height = self.size

            insert = (i*self.ticks_gap_width, 0) # insert should be the up-left most
            if self.mirror:
                insert = (i*self.ticks_gap_width, -height/2.0)
            size = (self.label_width, height/2.0) # size should be positive

            tick_label = str(i)
            tick = Tick(insert, size, tick_label, 
                    anchor_object=line, self_pos='middle',
                    mirror=self.mirror).render()
            self._svg.add(tick)
   
        self._svg.add(line)
        self._svg.add(label_text)
        return self._svg 


class Tick(svgobject.Container):

    def __init__(self, insert, size, label, stroke='black', 
            rotate = 0, font_size=10, font_family='Arial',
            mirror = False, tick_self_width=3,
            **extra):

        super(Tick, self).__init__(insert, size, rotate=rotate, **extra)

        self.label = label
        self.stroke = stroke 
        self.font_size = font_size
        self.font_family = font_family
        self.mirror = mirror
        self.tick_self_width = tick_self_width

    def render(self):
        insert = 0, 0
        size = 0, self.tick_self_width 
        related_insert=('up','middle')
        if self.mirror:
            related_insert=('bottom','middle')

        tick = svgobject.Line(insert, size, stroke=self.stroke, transform=self.transform,
                parent=self._svg, related_insert=related_insert).render()

        related_insert = ('bottom', 'middle')
        if self.mirror:
            related_insert = ('up', 'middle')

        label_text = svgobject.Text(self.label, insert, anchor_object=tick, related_insert=related_insert, 
                font_size=self.font_size, font_family=self.font_family).render()

        self._svg.add(tick)
        self._svg.add(label_text)

        return self._svg
