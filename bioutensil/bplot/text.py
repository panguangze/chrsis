import os
from io import BytesIO
import matplotlib
import matplotlib.pyplot as plt

import svgwrite
from bioutensil.bplot import svgobject
from bioutensil.bplot import svgbase

plt.rc('mathtext', fontset='cm')
matplotlib.use('agg')


def measure_text(text, font_size=10):
    '''
    measure text length for Arial,
    ref: https://bl.ocks.org/tophtucker/62f93a4658387bb61e4510c37e2e97cf
    '''
    widths = [0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0.2796875,0.2765625,0.3546875,0.5546875,0.5546875,0.8890625,0.665625,0.190625,0.3328125,0.3328125,0.3890625,0.5828125,0.2765625,0.3328125,0.2765625,0.2765625,0.5546875,0.5546875,0.5546875,0.5546875,0.5546875,0.5546875,0.5546875,0.5546875,0.5546875,0.5546875,0.2765625,0.2765625,0.5828125,0.5828125,0.5828125,0.5546875,1.0140625,0.66875,0.665625,0.721875,0.721875,0.665625,0.609375,0.7765625,0.721875,0.2765625,0.5,0.665625,0.5546875,0.8328125,0.721875,0.7765625,0.665625,0.7765625,0.721875,0.665625,0.609375,0.721875,0.665625,0.94375,0.665625,0.665625,0.609375,0.2765625,0.2765625,0.2765625,0.46875,0.58125,0.3328125,0.5546875,0.5546875,0.5,0.5546875,0.5546875,0.3125,0.5546875,0.5546875,0.221875,0.2671875,0.5,0.221875,0.8328125,0.5546875,0.5546875,0.5546875,0.5546875,0.3453125,0.5,0.2765625,0.5546875,0.5,0.721875,0.5,0.5,0.5,0.3328125,0.259375,0.3328125,0.5828125]
    avg = 0.5273190789473683
    temps = [widths[ord(c)] if ord(c) < len(widths) else avg for c in text]
    l = sum(temps) * font_size
    return l


def tex2svg(latex, fontsize=100, dpi=300, output=BytesIO()):
    fig = plt.figure(figsize=(0.01, 0.01))
    fig.text(0, 0, latex, fontsize=fontsize)
    fig.savefig(output, dpi=dpi, transparent=True, format='svg',
                bbox_inches='tight', pad_inches=0.01, frameon=False)
    plt.close(fig)
    if isinstance(output, BytesIO):
        output.seek(0)
        return output.read()
    else:
        return


class LatexText(svgbase.SVGObject):

    def __init__(self, text, name, insert=(0, 0), size=(5, 5),
                 font_size=10,
                 anchor_object=None, anchor_pos='',
                 target_object=None, target_pos='',
                 **kwargs):

        self.x, self.y = insert
        self.width, self.height = size

        super(LatexText, self).__init__(insert, size,
                                        font_size=font_size,
                                        anchor_object=anchor_object,
                                        anchor_pos=anchor_pos,
                                        target_object=target_object,
                                        target_pos=target_pos,
                                        **kwargs)

        if not os.path.exists('latex'):
            os.makedirs('latex')
        fn = 'latex/{}.svg'.format(name)
        tex2svg(text, fontsize=font_size, dpi=300,
                output=fn)

        self.svgobj = svgwrite.image.Image(fn, self.insert, self.size,
                                           **kwargs)


class TextArea(svgobject.Container):

    def __init__(self, text, insert,
                 font_size=10,
                 **kwargs):
        self.texts = text.split('\n')
        self.font_size = font_size

        size = self._get_size()
        super(TextArea, self).__init__(insert, size, **kwargs)
        self.split_items = sum(self.split(split_pattern=(1, len(self.texts))), [])

    def render(self):
        for text, text_container in zip(self.texts, self.split_items):
            text_container = text_container.render()
            insert = 0, 0
            text_object = svgobject.Text(text, insert,
                                         parent=text_container, related_insert='center:middle').render()

            text_container.add(text_object)
            self._svg.add(text_container)

        return self._svg

    def _get_size(self):
        width = max(len(x) for x in self.texts)*8
        height = len(self.texts)*self.font_size
        size = width, height
        return size
