import numpy as np

def generate_naca4_airfoil(m=2, p=4, t=12, n=128, alpha_deg=0):
    """
    Generate a rotated (angle of attack) smooth NACA 4-digit airfoil.
    Returns upper and lower surfaces as closed loop.
    """
    m /= 100
    p /= 10
    t /= 100
    alpha = np.radians(alpha_deg)

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

    x_full = np.concatenate([xu[::-1], xl[1:]])
    y_full = np.concatenate([yu[::-1], yl[1:]])

    # Rotate by angle of attack
    x_rot = x_full * np.cos(alpha) - y_full * np.sin(alpha)
    y_rot = x_full * np.sin(alpha) + y_full * np.cos(alpha)

    return x_rot, y_rot, x_full, y_full

def write_geo_file_with_three_lines(x, y, filename="naca_airfoil.geo"):
    lines = ["// GMSH .geo file for NACA airfoil with 3 lines\n"]
    for i, (xi, yi) in enumerate(zip(x, y)):
        lines.append(f"Point({i+1}) = {{{xi:.8f}, {yi:.8f}, 0.0, 1.0}};\n")

    # Find division points at x = 0.05 (5% chord from leading edge)
    x_arr = np.array(x)
    front_upper = np.argmax(x_arr <= 0.05)
    front_lower = len(x) - np.argmax(x_arr[::-1] <= 0.05) - 1

    # Define 3 line segments
    lines.append(f"Spline(1) = {{{1}:{front_upper+1}}};\n")
    lines.append(f"Spline(2) = {{{front_upper+1}:{front_lower+1}}};\n")
    lines.append(f"Spline(3) = {{{front_lower+1}:{len(x)}}};\n")

    with open(filename, 'w') as f:
        f.writelines(lines)
    print(f"Written GMSH .geo file with 3 lines to: {filename}")

if __name__ == '__main__':
    x, y, _, _ = generate_naca4_airfoil(m=2, p=4, t=12, n=128, alpha_deg=5)
    write_geo_file_with_three_lines(x, y, "naca2412_airfoil.geo")