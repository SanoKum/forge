import numpy as np

def generate_naca4_airfoil(m=2, p=4, t=12, n=255, alpha_deg=0):
    m /= 100
    p /= 10
    t /= 100
    alpha = np.radians(alpha_deg)

    beta = np.linspace(0, np.pi, n)
    x = (1 - np.cos(beta)) / 2

    yt = 5 * t * (0.2969*np.sqrt(x) - 0.1260*x - 0.3516*x**2 + 0.2843*x**3 - 0.1015*x**4)

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

    xu[0], yu[0] = (xu[0] + xl[0]) / 2, (yu[0] + yl[0]) / 2
    xl[0], yl[0] = xu[0], yu[0]

    x_full = np.concatenate([xu[::-1], xl[1:]])
    y_full = np.concatenate([yu[::-1], yl[1:]])

    x_rot = x_full * np.cos(alpha) - y_full * np.sin(alpha)
    y_rot = x_full * np.sin(alpha) + y_full * np.cos(alpha)

    return x_rot, y_rot

def write_geo_file_with_outer_box_and_mesh(x, y, filename="naca2412_with_box.geo"):
    lines = ["// GMSH .geo file with airfoil and outer box\n"]
    for i, (xi, yi) in enumerate(zip(x, y)):
        lines.append(f"Point({i+1}) = {{{xi:.8f}, {yi:.8f}, 0.0, 0.001}};\n")

    num_pts = len(x)
    spline_pts = ', '.join(str(i+1) for i in range(num_pts))
    lines.append(f"Spline(1) = {{{spline_pts}}};\n")
    lines.append(f"Line(2) = {{{num_pts}, 1}};\n")
    lines.append("Line Loop(10) = {1, 2};\n")

    # Outer box (10 chord lengths wide)
    box_start = num_pts + 1
    lines += [
        f"Point({box_start}) = {{-5, -5, 0, 0.5}};\n",
        f"Point({box_start+1}) = {{10, -5, 0, 0.5}};\n",
        f"Point({box_start+2}) = {{10, 5, 0, 0.5}};\n",
        f"Point({box_start+3}) = {{-5, 5, 0, 0.5}};\n",
        f"Line(20) = {{{box_start}, {box_start+1}}};\n",
        f"Line(21) = {{{box_start+1}, {box_start+2}}};\n",
        f"Line(22) = {{{box_start+2}, {box_start+3}}};\n",
        f"Line(23) = {{{box_start+3}, {box_start}}};\n",
        "Line Loop(30) = {20, 21, 22, 23};\n",
        "Plane Surface(1) = {10, 30};\n"
    ]

    lines.append("Mesh.CharacteristicLengthMin = 0.001;\n")
    lines.append("Mesh.CharacteristicLengthMax = 0.5;\n")

    lines += [
        "Field[1] = BoundaryLayer;\n",
        "Field[1].EdgesList = {1, 2};\n",
        "Field[1].hwall_n = 0.0005;\n",
        "Field[1].thickness = 0.01;\n",
        "Field[1].ratio = 1.2;\n",
        "Field[1].nLayers = 3;\n",
        "BoundaryLayer Field = 1;\n",
        "Extrude {0, 0, 0.01} {\n",
        "  Surface{1}; Layers {1}; Recombine;\n",
        "}\n",
        "Physical Surface(\"inlet\", 1) = {49};\n",
        "Physical Surface(\"outlet\", 2) = {57};\n",
        "Physical Surface(\"top\", 3) = {53};\n",
        "Physical Surface(\"bot\", 4) = {61};\n",
        "Physical Surface(\"wall\", 5) = {41};\n",
        "Physical Surface(\"side1\", 6) = {41};\n",
        "Physical Surface(\"side2\", 7) = {41};\n",
        "Physical Surface(\"fluid\", 8) = {1};\n"
    ]

    with open(filename, 'w') as f:
        f.writelines(lines)
    print(f"Written GMSH .geo file with outer box and boundary layer to: {filename}")

if __name__ == '__main__':
    x, y = generate_naca4_airfoil(m=2, p=4, t=12, n=255, alpha_deg=5)
    write_geo_file_with_outer_box_and_mesh(x, y, "naca2412_with_box.geo")