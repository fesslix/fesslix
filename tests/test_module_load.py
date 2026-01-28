#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Python script to test if the Fesslix module can be imported

Copyright (C) 2010-2026 Wolfgang Betz

Copying and distribution of this file, with or without modification,
are permitted in any medium without royalty provided the copyright
notice and this notice are preserved.  This file is offered as-is,
without any warranty.
"""

import fesslix as flx
print(flx.__version__)
flx.load_engine()

import fesslix.gpr
import fesslix.tools
import fesslix.model_templates


import numpy as np

flx.print_info()

