#pragma once

#include "IAlgorithm.h"
using IAlgorithm = ses::IAlgorithm;

IAlgorithm ses::create_cg_algorithm();
IAlgorithm ses::create_pcg_algorithm();

