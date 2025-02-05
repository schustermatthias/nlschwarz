import numpy as np


def load_function_1(x):
    return 5.0


def load_function_2(x):
    return 5.0


def load_function_3(x):
    return 1.0


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
        "function": "schwarz_dirichlet_1",
        "fractional_s": 0.5,
        "horizon": 0.1,
        "outputdim": 1
    },
    {
        "function": "schwarz_dirichlet_2",
        "fractional_s": 0.5,
        "horizon": 0.1,
        "outputdim": 1
    },
    {
        "function": "schwarz_dirichlet_5",
        "fractional_s": 0.5,
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
          'function':load_function_1},
         {'weights': dx,
          'points': Px,
          'function': load_function_2},
         {'weights': dx,
          'points': Px,
          'function': load_function_3}
         ]

