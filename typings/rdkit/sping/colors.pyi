from typing import Any

class Color:
    def __init__(self, red: int = ..., green: int = ..., blue: int = ...) -> None: ...
    def __setattr__(self, name: Any, value: Any) -> None: ...
    def __mul__(self, x: Any): ...
    def __rmul__(self, x: Any): ...
    def __truediv__(self, x: Any): ...
    def __div__(self, x: Any): ...
    def __rdiv__(self, x: Any): ...
    def __add__(self, x: Any): ...
    def __sub__(self, x: Any): ...
    def __hash__(self) -> Any: ...
    def __cmp__(self, other: Any): ...
    def toHexRGB(self): ...
    def toHexStr(self): ...

def HexColor(val: Any): ...

aliceblue: Any
antiquewhite: Any
aqua: Any
aquamarine: Any
azure: Any
beige: Any
bisque: Any
black: Any
blanchedalmond: Any
blue: Any
blueviolet: Any
brown: Any
burlywood: Any
cadetblue: Any
chartreuse: Any
chocolate: Any
coral: Any
cornflower: Any
cornsilk: Any
crimson: Any
cyan: Any
darkblue: Any
darkcyan: Any
darkgoldenrod: Any
darkgray: Any
darkgreen: Any
darkkhaki: Any
darkmagenta: Any
darkolivegreen: Any
darkorange: Any
darkorchid: Any
darkred: Any
darksalmon: Any
darkseagreen: Any
darkslateblue: Any
darkslategray: Any
darkturquoise: Any
darkviolet: Any
deeppink: Any
deepskyblue: Any
dimgray: Any
dodgerblue: Any
firebrick: Any
floralwhite: Any
forestgreen: Any
fuchsia: Any
gainsboro: Any
ghostwhite: Any
gold: Any
goldenrod: Any
gray: Any
grey = gray
green: Any
greenyellow: Any
honeydew: Any
hotpink: Any
indianred: Any
indigo: Any
ivory: Any
khaki: Any
lavender: Any
lavenderblush: Any
lawngreen: Any
lemonchiffon: Any
lightblue: Any
lightcoral: Any
lightcyan: Any
lightgoldenrodyellow: Any
lightgreen: Any
lightgrey: Any
lightpink: Any
lightsalmon: Any
lightseagreen: Any
lightskyblue: Any
lightslategray: Any
lightsteelblue: Any
lightyellow: Any
lime: Any
limegreen: Any
linen: Any
magenta: Any
maroon: Any
mediumaquamarine: Any
mediumblue: Any
mediumorchid: Any
mediumpurple: Any
mediumseagreen: Any
mediumslateblue: Any
mediumspringgreen: Any
mediumturquoise: Any
mediumvioletred: Any
midnightblue: Any
mintcream: Any
mistyrose: Any
moccasin: Any
navajowhite: Any
navy: Any
oldlace: Any
olive: Any
olivedrab: Any
orange: Any
orangered: Any
orchid: Any
palegoldenrod: Any
palegreen: Any
paleturquoise: Any
palevioletred: Any
papayawhip: Any
peachpuff: Any
peru: Any
pink: Any
plum: Any
powderblue: Any
purple: Any
red: Any
rosybrown: Any
royalblue: Any
saddlebrown: Any
salmon: Any
sandybrown: Any
seagreen: Any
seashell: Any
sienna: Any
silver: Any
skyblue: Any
slateblue: Any
slategray: Any
snow: Any
springgreen: Any
steelblue: Any
tan: Any
teal: Any
thistle: Any
tomato: Any
turquoise: Any
violet: Any
wheat: Any
white: Any
whitesmoke: Any
yellow: Any
yellowgreen: Any
transparent: Any
