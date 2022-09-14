from q2_types.per_sample_sequences import PairedEndSequencesWithQuality
from q2_types.sample_data import SampleData
from qiime2.core.type import Str
from qiime2.plugin import Plugin

from q2_pan_classifier.actions.prep_sequences import import_reads

from . import __version__

# This is the plugin object. It is what the framework will load and what an
# interface will interact with. Basically every registration we perform will
# involve this object in some way.

plugin = Plugin(name="pan-classifier",
                version=__version__,
                package="q2_pan_classifier",
                website="https://github.com/CRideTGen/q2-pan-classifier")

plugin.methods.register_function(
    function=import_reads,
    inputs={},
    parameters={"sequence_directory": Str},
    outputs=[('sequences_out', SampleData[PairedEndSequencesWithQuality])],
    input_descriptions={},
    parameter_descriptions={"sequence_directory": "Path to sequencing directory"},
    output_descriptions={'sequences_out': 'Bowtie2 index.'},
    name='Import read data using manifest',
    description='Import read data using manifest protocol.'
)