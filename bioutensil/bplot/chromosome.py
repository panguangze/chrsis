
import svgobject
import shape

class Chrom(svgobject.Container):

    def __init__(self, insert, size,
            chrom, p_start, p_end, q_start, q_end,
            chrom_width=3, rx=5, ry=5,
            font_size=10, font_family='Arial',
            h_margin=80, v_margin=5, cytoband=None,
            opacity=0.8,
            **extra):


        self.chrom = chrom
        self.p_start = p_start
        self.p_end = p_end
        self.q_start = q_start
        self.q_end = q_end

        self.cytoband_info=cytoband
        self.cytoband_color_scheme={
                'stalk': '#f0f0f0',
                'gvar': '#f0f0f0',
                'acen': '#960015', #'#f0f0f0',
                'gneg': '#d9d9d9',
                'gpos25': '#bdbdbd',
                'gpos50': '#969696',
                'gpos75': '#737373',
                'gpos100': '#525252',
        }

        self.chrom_width = chrom_width
        self.rx = rx
        self.ry = ry
        self.font_size = font_size
        self.font_family = font_family
        self.h_margin=h_margin
        self.v_margin=v_margin

        size = (q_end, chrom_width)
        super(Chrom, self).__init__(insert, size, h_margin=h_margin, v_margin=v_margin, **extra)


    def render(self):
        anchor = self._render_cytoband()
        self._render_label(anchor)
        return self._svg

    def _render_label(self, anchor):
        # label
        label_insert = (None, None)
        related_insert = ('center', 'left')
        label = svgobject.Text(self.chrom, label_insert, font_size=self.font_size, font_family=self.font_family, anchor_object=anchor, related_insert=related_insert).render()
        self._svg.add(label)

    def _render_arm(self):
        # p_arm
        p_insert = (self.x + self.h_margin + self.p_start, 0)
        p_size = (self.p_end - self.p_start, self.chrom_width)
        q_insert = (self.x + self.h_margin + self.q_start, 0)
        q_size = (self.q_end - self.q_start, self.chrom_width)
        for arm_name, insert, size, arm_start in [('q', q_insert, q_size, self.q_start), ('p', p_insert, p_size, self.p_start)]:
            arm = Arm(insert, size, self.chrom, arm_name, arm_start,
                    cytoband=self.cytoband_info,
                    rx=self.rx, ry=self.ry,
                    parent=self._svg, related_insert='center').render()
            self._svg.add(arm)
        return arm

    def _render_cytoband(self):
        insert = self.x + self.h_margin, 0
        size = self.q_end - self.p_start, self.chrom_width
        cytoband = CytoBand(insert, size, self.chrom,
                opacity=self.opacity,
                cytoband=self.cytoband_info, cytoband_color_scheme=self.cytoband_color_scheme,
                parent=self._svg, related_insert='center').render()
        self._svg.add(cytoband)

        return cytoband

class CytoBand(svgobject.Container):
    def __init__(self, insert, size,
            chrom, cytoband=None, cytoband_color_scheme=None,
            **extra):

        super(CytoBand, self).__init__(insert, size, **extra)
        self.chrom = chrom

        self.cytoband_info=cytoband
        if not cytoband_color_scheme:
            cytoband_color_scheme= {
                'stalk': '#f0f0f0',
                'gvar': '#f0f0f0',
                'acen': '#f0f0f0',
                'gneg': '#d9d9d9',
                'gpos25': '#bdbdbd',
                'gpos50': '#969696',
                'gpos75': '#737373',
                'gpos100': '#525252',
            }
        self.cytoband_color_scheme = cytoband_color_scheme



    def render(self):
        width, height = self.size
        cytoband_info = self.cytoband_info['cytoband']
        chrom2endposition = self.cytoband_info['chrom2endposition']
        for chrom, position, start, end, color  in cytoband_info:
            if chrom != self.chrom: continue
            if color == 'acen':
                if 'p' in position:
                    points = [(start, 0), (start, height), (end, 2.0*height/3), (end, 1.0*height/3)]
                if 'q' in position:
                    points = [(start, 1.0*height/3), (start, 2.0*height/3), (end, height), (end, 0)]
                stroke_pattern = [0,1,0,1]
                cytoband = shape.Polygon(points,
                        stroke_width=0.2, stroke_pattern=stroke_pattern,
                        opacity=self.opacity,
                        fill=self.cytoband_color_scheme[color],
                        ).render()
            else:
                insert = start, 0
                size = end - start, height
                front, end = chrom2endposition[chrom]
                if position == front:   stroke_lines = ['up','bottom','left']
                elif position == end:   stroke_lines = ['up','bottom','right']
                else:                   stroke_lines = ['up','bottom']
                cytoband = svgobject.Rect(insert, size,
                        stroke_width=0.3, stroke_lines=stroke_lines,
                        opacity=self.opacity,
                        fill=self.cytoband_color_scheme[color],
                        ).render()
            self._svg.add(cytoband)
        return self._svg


class Arm(svgobject.Container):
    def __init__(self, insert, size,
            chrom, name='p', arm_start=0,
            cytoband=None,
            rx=5, ry=5, opacity=0.5,
            **extra):

        super(Arm, self).__init__(insert, size, **extra)
        self.chrom = chrom
        self.name = name
        self.arm_start = arm_start

        self.cytoband_info=cytoband
        self.cytoband_color_scheme={
                'stalk': '#f0f0f0',
                'gvar': '#f0f0f0',
                'acen': '#f0f0f0',
                'gneg': '#d9d9d9',
                'gpos25': '#bdbdbd',
                'gpos50': '#969696',
                'gpos75': '#737373',
                'gpos100': '#525252',
        }

        self.rx = rx
        self.ry = ry
        self.opacity = opacity


    def render(self):
        self._render_cytoband()
        insert = (0,0)
        self.arm = svgobject.Rect(insert, self.size,
                rx=self.rx, ry=self.ry, opacity=self.opacity).render()
        self._svg.add(self.arm)
        return self._svg

    def _render_cytoband(self):
        width, height = self.size
        for chrom, position, start, end, color  in self.cytoband_info:
            if chrom != self.chrom: continue
            if self.name not in position: continue
            insert = start - self.arm_start , 0
            size = end - start, height
            rx = 0
            ry = 0
            cytoband = svgobject.Rect(insert, size, rx=rx, ry=ry,
                    fill=self.cytoband_color_scheme[color], stroke='white',
                    parent=self._svg, related_insert='center').render()
            self._svg.add(cytoband)



