"""メッシュ層: 構造化 2D quad 生成 → gmsh msh4.1 直書き → convertGmshToForge。"""

from .mesh2d import Mesh2DParams, generate_axisym_mesh, write_msh41_2d

__all__ = ["Mesh2DParams", "generate_axisym_mesh", "write_msh41_2d"]
