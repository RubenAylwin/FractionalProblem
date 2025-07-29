//////////////////////////////
// Unit tests for the code. //
//////////////////////////////

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE test suite
#include <boost/test/included/unit_test.hpp>
#include <Msg.h>
useMessages("TEST_DEBUG");
#include "geomTest.h"
#include "discreteTest.h"
#include "scalarFunctionTest.h"
#include "integrationTest.h"
#include "greenTest.h"
#include "utilsTest.h"
#include "empiricalInterpolationTest.h"
#include "quasiMonteCarloTest.h"
#include "riemannLiouvilleTest.h"
#include "meshTest.h"
#include "discreteMeshTest.h"
