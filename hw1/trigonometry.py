import sys
import numpy as np
import matplotlib.pyplot as plt

# Note: I am trying to be more efficient in my coding since its comp physics and we kind need efficiency
# even if realistically I would be doing C++ if I want efficiency. Just good practice I suppose.
# More accurately what Im practicing is not exactly efficiency and more numerical stability but both are good ig.

VALID_FLAGS = {"--function", "--write", "--read_from_file", "--print"}
VALID_FMTS = {"jpeg", "jpg", "eps", "pdf", "png"}

def sinc(x, out = None):
    """
    I tried to write sinc in a optimized manner. 
    I further bound the error to below |x| < 1E-4 via the second order taylore expansion.
    This truncation error is bound by x⁴/120 ≈ 8.33E-19, which is well below the machine delta of a double (≈2.22E-16).
    Therefore, there will be no indeterminate form nor floating point errors near 0
    """
    x = np.asarray(x, dtype = np.float64)
    if out is None:
        out = np.empty_like(x)
        
    near_zero = np.abs(x) < 1e-4
    far = ~near_zero # Here, we used the bitwise inversion 
    
    x_near = x[near_zero]
    out[near_zero] = 1.0 - (x_near * x_near) / 6 # This is the taylort expansion
    
    x_far = x[far]
    out[far] = np.sin(x_far) / x_far
        
    return out

VALID_FUNCS = {"cos": np.cos, "sin": np.sin, "sinc": sinc}

def parse_args():
    # Here, we process our arguments into a large list
    tokens = []
    for arg in sys.argv[1:]:
        if arg.startswith("--") and "=" in arg:
            flag, val = arg.split("=",1)
            tokens.extend([flag, val])
        else:
            tokens.append(arg)
    
    # Here, we will map our tokens to flags
    opts = {}
    current_flag = None
    for token in tokens:
        if token.startswith("--"):
            if token not in VALID_FLAGS:
                sys.exit(f"Error: Unknown flag '{token}'")
            current_flag = token
            opts.setdefault(current_flag, [])
        elif current_flag:
            opts[current_flag].extend(token.replace(",", " ").split())
        else:
            sys.exit(f"Error: Argument '{token}'lacks an associated flag")
            
    # We must now make sure all the minimally required flags exist and also check all existing flags have their output:
    if "--read_from_file" in opts:
        if not opts["--read_from_file"]:
            sys.exit("Error: '--read_from_file' requires a filename.")
    elif "--function" in opts:
        if not opts["--function"]:
            sys.exit("Error: '--function' requires at least one value")
        for f in opts["--function"]:
            if f not in VALID_FUNCS:
                sys.exit(f"Error: '{f}' is not a supported function. Choose from: 'sin', 'cos', 'sinc'.")
        opts["--function"] = list(dict.fromkeys(opts["--function"])) 
        # While we could do list(set(...)) above to remove duplicates, this would not preserve plotting order.
    else:
        sys.exit("Error: Must provide either '--read_from_file' or '--function'")

    if "--write" in opts and not opts["--write"]:
        sys.exit("Error: '--write' requires a output filename.")
        
    if "--print" in opts:
        if not opts["--print"]:
            sys.exit("Error: '--print' requires at least one format.")
        for fmt in opts["--print"]:
            if fmt.lower() not in VALID_FMTS:
                sys.exit(f"Error: '{fmt}' is not a supported format.")
        
    return opts


def main():
    opts = parse_args()

    # Here do we parts (a) and (c)
    if "--read_from_file" in opts:
        filename = opts["--read_from_file"][0]
        try:
            with open(filename, "r") as f:
                # I do believe some NumPy stuff causes issues with a '#' so Im just safely splitting them since idk how to handle them
                headers = f.readline().replace("#", "").split()[1:]
            data = np.loadtxt(filename, skiprows=1, ndmin=2)
            x = data[:, 0]
            curves = {name: data[:, i + 1] for i, name in enumerate(headers)}
        except Exception as err:
            sys.exit(f"Error reading '{filename}': {err}")
    else:
        x = np.linspace(-10, 10, 401)
        curves = {name: VALID_FUNCS[name](x) for name in opts["--function"]}

    # Here be do part (b)
    if "--write" in opts:
        filename = opts["--write"][0]
        header = "x " + " ".join(curves.keys())
        matrix = np.column_stack([x, *curves.values()])
        np.savetxt(filename, matrix, header=header, comments="", fmt="%16.8e")

    # This is obviously the plotting:
    plt.figure(figsize=(8, 4.5))
    for name, y in curves.items():
        plt.plot(x, y, label=f"{name}(x)", linewidth=1.5)

    plt.title("Plot of Trigonometric Functions")
    plt.xlabel("x")
    plt.ylabel("f(x)")
    plt.axhline(0, color="black", linestyle="--", linewidth=0.6, alpha=0.5)
    plt.axvline(0, color="black", linestyle="--", linewidth=0.6, alpha=0.5)
    plt.grid(True, linestyle=":", alpha=0.6)
    plt.legend()
    plt.tight_layout()

    # Here we do part (d)
    if "--print" in opts:
        for fmt in opts["--print"]:
            plt.savefig(f"plot.{fmt.lower()}", format=fmt.lower(), dpi=300)
    else:
        plt.show()


if __name__ == "__main__":
    main()