#!/usr/bin/env python
import tandem
seq = "AAGAAGAAGATGAAGAGAAGTTTTT"
print tandem.search_issr(seq, 3, 8, 3, 1, 2, 10, 50)
