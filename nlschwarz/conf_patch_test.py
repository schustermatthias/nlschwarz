import numpy as np

def load_function(x):
    # 2
    return 0.0


# Px = np.array([[0.33333333, 0.33333333],
#               [0.45929259, 0.45929259],
#               [0.45929259, 0.08141482],
#               [0.08141482, 0.45929259],
#               [0.17056931, 0.17056931],
#               [0.17056931, 0.65886138],
#               [0.65886138, 0.17056931],
#               [0.05054723, 0.05054723],
#               [0.05054723, 0.89890554],
#               [0.89890554, 0.05054723],
#               [0.26311283, 0.72849239],
#               [0.72849239, 0.00839478],
#               [0.00839478, 0.26311283],
#               [0.72849239, 0.26311283],
#               [0.26311283, 0.00839478],
#               [0.00839478, 0.72849239]])

# dx = 0.5 * np.array([0.14431560767779,
#                      0.09509163426728,
#                      0.09509163426728,
#                      0.09509163426728,
#                      0.10321737053472,
#                      0.10321737053472,
#                      0.10321737053472,
#                      0.03245849762320,
#                      0.03245849762320,
#                      0.03245849762320,
#                      0.02723031417443,
#                      0.02723031417443,
#                      0.02723031417443,
#                      0.02723031417443,
#                      0.02723031417443,
#                      0.02723031417443])

Px = np.array([[0.33333333333333,    0.33333333333333],
                  [0.47014206410511,    0.47014206410511],
                  [0.47014206410511,    0.05971587178977],
                  [0.05971587178977,    0.47014206410511],
                  [0.10128650732346,    0.10128650732346],
                  [0.10128650732346,    0.79742698535309],
                  [0.79742698535309,    0.10128650732346]])
dx = 0.5 * np.array([0.22500000000000,
                     0.13239415278851,
                     0.13239415278851,
                     0.13239415278851,
                     0.12593918054483,
                     0.12593918054483,
                     0.12593918054483])

Py = Px  # np.array([[0.33333333, 0.33333333]])
dy = dx  # 0.5 * np.array([1.0])

kernels = [{
        "function": "constant",
        "fractional_s": 0.6,
        "horizon": 0.1,
        "outputdim": 1
    },
    {
        "function": "fractional",
        "fractional_s": 0.6,
        "horizon": 0.1,
        "outputdim": 1
    }
]

confs = [{
        # "savePath": "pathA",
        "ansatz": "CG", #DG
        "is_fullConnectedComponentSearch": 0,
        "approxBalls": {
            "method": "retriangulate",
            "isPlacePointOnCap": True,  # required for "retriangulate" only
            #"averageBallWeights": [1., 1., 1.]  # required for "averageBall" only
        },
        "closeElements": "fractional",
        "quadrature": {
            "outer": {
                "points": Px,
                "weights": dx
            },
            "inner": {
                "points": Px,
                "weights": dx
            },
            "tensorGaussDegree": 5,  # Degree of tensor Gauss quadrature for weakly singular kernels.
        },
        "verbose": False
    },
{
        # "savePath": "pathA",
        "ansatz": "CG", #DG
        "is_fullConnectedComponentSearch": 0,
        "approxBalls": {
            "method": "retriangulate",
            "isPlacePointOnCap": True,  # required for "retriangulate" only
            #"averageBallWeights": [1., 1., 1.]  # required for "averageBall" only
        },
        "closeElements": "fractional",
        "quadrature": {
            "outer": {
                "points": Px,
                "weights": dx
            },
            "inner": {
                "points": Px,
                "weights": dx
            },
            "tensorGaussDegree": 5,  # Degree of tensor Gauss quadrature for weakly singular kernels.
        },
        "verbose": False
    }
]

loads = [{'weights': dx,
          'points': Px,
          'function':load_function},
         {'weights': dx,
          'points': Px,
          'function': load_function}
         ]
