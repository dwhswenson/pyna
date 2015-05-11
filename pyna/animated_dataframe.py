# parts stolen from: http://jakevdp.github.io/blog/2013/05/12/embedding-matplotlib-animations/

import matplotlib.animation as anim
import matplotlib.pyplot as plt
import numpy as np

from tempfile import NamedTemporaryFile

VIDEO_TAG = """<video controls>
 <source src="data:video/x-m4v;base64,{0}" type="video/mp4">
 Your browser does not support the video tag.
</video>"""

def anim_to_html(anim):
    if not hasattr(anim, '_encoded_video'):
        with NamedTemporaryFile(suffix='.mp4') as f:
            anim.save(f.name, fps=20, extra_args=['-vcodec', 'libx264'])
            video = open(f.name, "rb").read()
        anim._encoded_video = video.encode("base64")
    
    return VIDEO_TAG.format(anim._encoded_video)


from IPython.display import HTML

def display_animation(anim):
    plt.close(anim._fig)
    return HTML(anim_to_html(anim))


class AnimatedDataFrame(object):
    def __init__(self, df, columns=None):
        self.df = df
        if columns != None:
            self.xvals = columns
        else:
            self.xvals = df.columns.values
        self.fig, self.ax = plt.subplots()
        self.line, = self.ax.plot(self.xvals, self.df.loc[0,:])
        init = lambda : self.init()
        animate = lambda i : self.animate(i)
        self.ani = anim.FuncAnimation(self.fig, animate, self.df.index.values, init_func=init, interval=25, blit=True)

        
    def animate(self, i):
        self.line.set_ydata(self.df.loc[i,:])
        return self.line,
        
    def init(self):
        self.line.set_ydata(np.ma.array(self.xvals, mask=True))
        return self.line,
    
    def show(self):
        return display_animation(self.ani)

