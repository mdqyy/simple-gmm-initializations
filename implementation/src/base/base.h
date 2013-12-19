#ifndef BASE_H
#define BASE_H

#define NDEBUG

#include <cassert>
#include <cstddef>
#include <Eigen/Dense>
#include <random>
#include <stdint.h>

typedef double fp_type;
typedef Eigen::MatrixXd Matrix;
typedef Eigen::VectorXd Vector;
typedef Eigen::DenseIndex idx_type;

enum TestLabMode { DEFAULT_TEST, GENERATE_CSV };

enum GenMethod
{
	UNIFORM_MEANS_AND_RESTRICTED_COVAR
};

enum MeasurementError{ NO_ME, UNIFORM_NOISE };

enum AlgorithmID 
{ 
  // Learning Gaussian Mixture Models
  EM
};
  
enum InitMethod 
{ 
  EMPTY, 
  
  // GMM-intializations
  // ... dasgupta&schulman variants
  UNIFORM_SPHERICAL, 
  UNIFORM_SPHERICAL_WITH_PRUNING, 
  ADAPTIVE_SPHERICAL, 
  EXPERIMENTAL_ADAPTIVE, 
  // ... kmeans initializations "2gmm"
  UNIFORM_MEANS2GMM, 
  ADAPTIVE_MEANS2GMM, 
  GONZALEZ2GMM,
  // ... sequential algorithms
  ALTERNATELY_ADAPTIVEMEANS_MEANS2GMM, 
  GONZALEZ_FORGMM,
  GONZALEZ_KWEDLO,
  // ... restart_init
  RESTART_INIT,
  // .. agglomerative
  SAMPLE_AVGLINK_MEANS2GMM,   
};

enum AntiPropMode { INVERSEDENSITY, MAHALANOBIS };
enum PropMode { DENSITY, NEIGHBORHOOD };
enum CostMeasure{ NLL};


#ifndef M_PI
#define M_PI           3.14159265358979323846
#endif

#define gmmlab_throw(message) \
{printf("ERROR: %s\n     in file %s at line %d\n", message,__FILE__,__LINE__); throw(1);}

#endif
