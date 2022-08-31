#   Copyright 2021 Chase Ridenour

#   Licensed under the Apache License, Version 2.0 (the "License");
#   you may not use this file except in compliance with the License.
#   You may obtain a copy of the License at
#
#       http://www.apache.org/licenses/LICENSE-2.0
#
#   Unless required by applicable law or agreed to in writing, software
#   distributed under the License is distributed on an "AS IS" BASIS,
#   WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
#   See the License for the specific language governing permissions and
#   limitations under the License.

from qiime2.plugin import SemanticType

from ._format import DNAFastaNCBIDirFormat, NCBIAccFileDirectoryFormat
from q2_pan_classifier.plugin_setup import plugin

DNAFastaNCBI = SemanticType('DNAFastaNCBI')
NCBIAccFile = SemanticType("NCBIAccFile")



plugin.register_semantic_types(DNAFastaNCBI, NCBIAccFile)
plugin.register_semantic_type_to_format(DNAFastaNCBI, DNAFastaNCBIDirFormat)
plugin.register_semantic_type_to_format(NCBIAccFile, NCBIAccFileDirectoryFormat)