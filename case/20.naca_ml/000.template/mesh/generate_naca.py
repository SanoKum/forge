import numpy as np

def generate_naca4_airfoil(m=2, p=4, t=12, n=128):
    """
    Generate a smooth NACA 4-digit airfoil (e.g., NACA2412) with cosine-spaced x-coordinates.
    Returns the concatenated upper and lower surfaces as a closed loop.
    """
    m /= 100
    p /= 10
    t /= 100

    # Cosine spacing for better resolution near LE and TE
    beta = np.linspace(0, np.pi, n)
    x = (1 - np.cos(beta)) / 2

    yt = 5 * t * (0.2969*np.sqrt(x) - 0.1260*x - 0.3516*x**2 
                  + 0.2843*x**3 - 0.1015*x**4)

    yc = np.where(x < p,
                  m / p**2 * (2*p*x - x**2),
                  m / (1-p)**2 * ((1 - 2*p) + 2*p*x - x**2))
    dyc_dx = np.where(x < p,
                      2*m / p**2 * (p - x),
                      2*m / (1-p)**2 * (p - x))
    theta = np.arctan(dyc_dx)

    xu = x - yt * np.sin(theta)
    yu = yc + yt * np.cos(theta)
    xl = x + yt * np.sin(theta)
    yl = yc - yt * np.cos(theta)

    # Remove duplicate trailing edge point
    x_coords = np.concatenate([xu[::-1], xl[1:]])
    y_coords = np.concatenate([yu[::-1], yl[1:]])

    return x_coords, y_coords

def write_geo_file(x, y, filename="naca_airfoil.geo"):
    lines = ["// GMSH .geo file for NACA airfoil\n"]
    for i, (xi, yi) in enumerate(zip(x, y)):
        lines.append(f"Point({i+1}) = {{{xi:.8f}, {yi:.8f}, 0.0, 1.0}};\n")
    for i in range(1, len(x)):
        lines.append(f"Line({i}) = {{{i}, {i+1}}};\n")
    lines.append(f"Line({len(x)}) = {{{len(x)}, 1}};\n")  # Close the loop

    with open(filename, 'w') as f:
        f.writelines(lines)
    print(f"Written GMSH .geo file to: {filename}")

if __name__ == '__main__':
    x, y = generate_naca4_airfoil(m=2, p=4, t=12, n=128)
    write_geo_file(x, y, "naca2412_airfoil.geo")
