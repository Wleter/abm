import signal
import sympy.physics.quantum.cg as cg
from sympy import S

def clebsch_gordan(double_j1: int, double_m1: int, double_j2: int, double_m2: int, double_j: int, double_m: int) -> float:
    signal.signal(signal.SIGINT, signal.SIG_DFL)
    result = cg.CG(S(double_j1)/2, S(double_m1)/2, S(double_j2)/2, S(double_m2)/2, S(double_j)/2, S(double_m)/2).doit()

    return float(result)

