"""Small module to load a colormap that's similar to tableau.

The main package (tableau20) takes no inputs and simply returns a nice
colormap that (at least from reading many blogs) is visually appealing to
a large number of people (and from personal experience is much better than
matplotlib defaults).

Colormap was adopted from a blog post on randalolson.com:
http://tinyurl.com/lum4ztm


"""


def tableau20(raw=None):
    """Return a vector of colors that mimics tableau."""
    # Create the tableau colors as RGB:
    tableau20 = [(31, 119, 180), (174, 199, 232), (255, 127, 14),
                 (255, 187, 120), (44, 160, 44), (152, 223, 138),
                 (214, 39, 40), (255, 152, 150),
                 (148, 103, 189), (197, 176, 213), (140, 86, 75),
                 (196, 156, 148), (227, 119, 194), (247, 182, 210),
                 (127, 127, 127), (199, 199, 199),
                 (188, 189, 34), (219, 219, 141), (23, 190, 207),
                 (158, 218, 229)]

    if raw is not None:
         return tableau20

    # Matplotlib assumes values that are scaled to 1
    for i in range(len(tableau20)):
        r, g, b = tableau20[i]
        tableau20[i] = (r / 255., g / 255., b / 255.)

    return tableau20
