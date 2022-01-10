#   Copyright 2021 Chase Ridenour
#
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

from q2_pan_classifier.actions.create_classifier import (create_classifier, generate_taxonomy, cross_validate_classifier)

from q2_pan_classifier.actions.prep_sequences import (prep_sequence_reads_single, prep_sequence_reads_paired)
from q2_pan_classifier.actions.classify_reads import (classify_reads, classify_reads_single, visualization_final)

# Good practice is to explicitely mark what is "available" for
# subpackage export. This is only used by the `from x import *` syntax, but
# it also keeps linters from complaining

__all__ = ['create_classifier', 'prep_sequence_reads_single','prep_sequence_reads_paired', 'classify_reads','classify_reads_single', 'visualization_final', 'generate_taxonomy', 'cross_validate_classifier']