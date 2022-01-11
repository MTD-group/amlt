
import numpy as np

def compute_vm_strain(B):
    return np.sqrt(2/3*np.linalg.norm(B)**2)


def random_sheer_matrix( target_von_mises_strain = 0.1, rng=None):
    if rng is None:
        rng = np.random.default_rng()

    # B is strain, M is the associated affine transformation matrix M = eye(3)+B

    ######
    # this is very ad hoc, but it makes the numbers look
    # nicer even before all the corrections
    scale = target_von_mises_strain/2.0
    B = np.zeros((3,3))
    for i in range(3):
        # subtracting trace()/3 from the diagonals changes their sdtev to 4/3 the stdev
        # before, this corrects that so all B_prime matrix elements have
        B [i,i] = rng.normal(scale = 3/4.* scale)
        for j in range(i):
            rv = rng.normal(scale = scale)
            B[i,j] = rv
            B[j,i] = rv

    #for i in range(3):
    #    B [i,i] = rng.normal(scale = 3/4.* scale)

    # Removing the trace for pure deviatoric strain
    B_prime = B - np.trace(B)/3 * np.eye(3)

    # Normalizing to the target von Mises strain
    vm_strain = compute_vm_strain(B_prime)
    B_prime_prime = (target_von_mises_strain/vm_strain ) * B_prime

    # removing any remainging traces of volume dilation
    M_prime_prime = B_prime_prime+np.eye(3)
    M_prime_prime_prime = M_prime_prime/(np.linalg.det(M_prime_prime))**(1/3)
    B_prime_prime_prime = M_prime_prime_prime - np.eye(3)

    return M_prime_prime_prime

