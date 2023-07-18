#pragma once
#include "IPreconditioner.h"

using IPreconditioner = ses::IPreconditioner;

IPreconditioner ses::create_sos_preconditioner();
