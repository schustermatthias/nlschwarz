import numpy as np

def load_function_1(x):
    return 10.0


def load_function_2(x):
    return 10.0


Px = np.array([[0.33333333, 0.33333333],
              [0.45929259, 0.45929259],
              [0.45929259, 0.08141482],
              [0.08141482, 0.45929259],
              [0.17056931, 0.17056931],
              [0.17056931, 0.65886138],
              [0.65886138, 0.17056931],
              [0.05054723, 0.05054723],
              [0.05054723, 0.89890554],
              [0.89890554, 0.05054723],
              [0.26311283, 0.72849239],
              [0.72849239, 0.00839478],
              [0.00839478, 0.26311283],
              [0.72849239, 0.26311283],
              [0.26311283, 0.00839478],
              [0.00839478, 0.72849239]])

dx = 0.5 * np.array([0.14431560767779,
                     0.09509163426728,
                     0.09509163426728,
                     0.09509163426728,
                     0.10321737053472,
                     0.10321737053472,
                     0.10321737053472,
                     0.03245849762320,
                     0.03245849762320,
                     0.03245849762320,
                     0.02723031417443,
                     0.02723031417443,
                     0.02723031417443,
                     0.02723031417443,
                     0.02723031417443,
                     0.02723031417443])


Py = np.array([[0.33333333, 0.33333333]])
dy = 0.5 * np.array([1.0])

kernels = [{
        "function": "constant",
        "fractional_s": 0.6,
        "horizon": 0.1,
        "outputdim": 1
    },
    {
        "function": "fractional",
        "fractional_s": 0.8,
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
          'function': load_function_1},
         {'weights': dx,
          'points': Px,
          'function': load_function_2}
         ]
