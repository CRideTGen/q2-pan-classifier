from qiime2.plugin import Plugin

from . import __version__

# This is the plugin object. It is what the framework will load and what an
# interface will interact with. Basically every registration we perform will
# involve this object in some way.

plugin = Plugin(name="pan-classifier",
                version=__version__,
                package="q2_pan_classifier",
                website="https://github.com/CRideTGen/q2-pan-classifier")
