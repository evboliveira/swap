import numpy as np
import glob, re, os

def find_folders(root):
    """Map N -> path of its NMentropy_swap.txt file, from subfolders named N<N>."""
    folders = {}
    for path in glob.glob(os.path.join(root, "N*")):
        if not os.path.isdir(path):
            continue
        m = re.fullmatch(r"N(\d+)", os.path.basename(path))
        if not m:
            continue
        N = int(m.group(1))
        fpath = os.path.join(path, "NMentropy_swap.txt")
        if os.path.isfile(fpath):
            folders[N] = fpath
    return dict(sorted(folders.items()))

def load_all(root):
    files = find_folders(root)
    if not files:
        raise FileNotFoundError(f"No N<N>/NMentropy_swap.txt files found under {root}")
    
    N_values = np.array(list(files.keys()))
    g_ref = None
    entropy = None
    error = None
    
    for i, N in enumerate(N_values):
        arr = np.loadtxt(files[N], comments='#')
        arr = arr[arr[:,0].argsort()]
        g, S, err = arr[:,0], arr[:,1], arr[:,2]
        
        if g_ref is None:
            g_ref = g
            entropy = np.full((len(N_values), len(g_ref)), np.nan)
            error   = np.full((len(N_values), len(g_ref)), np.nan)
        else:
            assert np.allclose(g, g_ref), f"g-grid mismatch in N={N}"
        
        entropy[i, :] = S
        error[i, :]   = err
    
    return N_values, g_ref, entropy, error

N_values, g_values, entropy, error = load_all('"results/28may2026/')

def f(N, g, tol=1e-6):
    """Return [entropy, error] for given N and g."""
    i = np.searchsorted(N_values, N)
    assert N_values[i] == N, f"N={N} not in dataset"
    j = np.argmin(np.abs(g_values - g))
    assert abs(g_values[j] - g) < tol, f"g={g} not on grid (closest {g_values[j]})"
    return np.array([entropy[i, j], error[i, j]])