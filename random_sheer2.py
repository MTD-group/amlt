import numpy as np


def deviatoric_strain_from_transform(A):
    strain=1/2*(A+A.T)-np.eye(3)
    dev = strain - np.eye(3)*np.trace(strain)/3
    return dev

def von_mises_from_transform(A):
    dev = deviatoric_strain_from_transform(A)
    vmeq = np.sqrt(2/3*(dev**2).sum() )
    return vmeq

def von_mises_of_sigma(sigma):
    delta_sigma = sigma-np.mean(sigma)
    vm = np.sqrt(2/3* np.sum(delta_sigma**2))
    return vm

def correct_sigma_product(sigma):
    ls = np.log(sigma)
    ls -= np.mean(ls)
    return np.exp(ls)

def correct_sigma_von_mises_for_target(sigma, vm_target):
    mean_sigma = np.mean(sigma)
    delta_sigma = sigma-mean_sigma
    correction_factor = vm_target/von_mises_of_sigma(sigma)
    return correction_factor*delta_sigma + mean_sigma


def random_sheer_matrix_svd(
    von_mises_strain,
    rng=None,
    verbose=False, atol=1e-12, imax=20):
    '''Create a random volume conserving affine transformation matrix with only 
    sheer components of specified von_mises_strain. '''

    vm_target = von_mises_strain
    if rng is None:
        random3x3 = np.random.rand(3,3)
    else:
        random3x3 = rng.random((3,3))

    randmat = (random3x3-0.5)*vm_target*2 + np.eye(3)

    W, S, Vh = np.linalg.svd(randmat, full_matrices=False)
    S0 = S.copy()
    # product(S)=1 >>> volume conserving transform
    # deviatoric is easy to define in principle axes (S) since there are no off diagonals
    # so how to satisfy product(S)=1 and von_mises=vm_target?

    # undoing the svd would be
    #transmat = W @ np.diagflat(S) @ Vh
    #but using only the polar part of the polar decomp
    #transmat = Vh.T @ np.diagflat(S) @ Vh
    vmeq = von_mises_of_sigma(S)
    i=0

    if verbose:
        unary_part = W @ Vh
        transmat = Vh.T @ np.diagflat(S) @ Vh
        print('vmeq (randmat)',von_mises_from_transform(randmat))
        print('det (randmat)', np.linalg.det(randmat))
        print('unary part',unary_part)
        print('polar part', transmat)
        print('vmeq of svd sigma', vmeq )
        print('vmeq', von_mises_from_transform(transmat))
        print('prod(S)', np.prod(S))
        print('det(transform)',np.linalg.det(transmat))
        print('--------corrections--------')
        print('step', i, S, '(vmeq) and error (%f) %f'%(vmeq, vmeq-vm_target))

    while np.abs(vmeq-vm_target) > atol and i<=imax:
        i+=1
        S=correct_sigma_von_mises_for_target(S,vm_target)
        S=correct_sigma_product(S)
        vmeq = von_mises_of_sigma(S)
        if verbose:
            print('step', i+1, S, '(vmeq) and error (%f) %e'%(vmeq, vmeq-vm_target))
    transmat = Vh.T @ np.diagflat(S) @ Vh

    if verbose:
        print('final matrix\n',transmat)
        print('vmeq', von_mises_from_transform(transmat))
        W, S, Vh = np.linalg.svd(transmat, full_matrices=False)
        polar_part = Vh.T @ np.diagflat(S) @ Vh
        unary_part = W @ Vh
        print('Polar decomp of final (sanity check)')
        print('unary part\n',unary_part)
        print('polar part\n', polar_part)

    return transmat


if __name__=='__main__':
    seed = 888
    rng = np.random.default_rng(seed)

    from ase.build import bulk
    from ase.visualize import view
    a=3.35 # not the real lattice parameter, but in ase gui it makes them touch slightly
    atoms0 = bulk('Cu',  'fcc', a=a, cubic=True).repeat((3,2,1))
    atoms0.translate(3*[a/4])
    atoms = atoms0.copy()


    sheermat = random_sheer_matrix_svd(von_mises_strain=0.05, rng=rng, verbose=True)
    new_cell = atoms.get_cell() @ sheermat
    atoms.set_cell(new_cell,scale_atoms=True)

    # from ase import io
    # io.write('before_and_after.traj',[atoms0,atoms])

    view([atoms0,atoms])

