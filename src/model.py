"""
doublet_core.py

Core data structures and functions for
'Automation of synthesis and ranking of cemented and air-spaced doublets'
(ITMO University, Romanova / Bakholdin, ~2022).

设计思路：
- 几乎每一个函数都对应流程图中的一个方块或一小段逻辑；
- 尽量避免深层 if 嵌套，采用“早返回 + 小函数”；
- 所有公式相关的部分都用 TODO 标出，后续按论文逐步补全。
"""

from __future__ import annotations

from dataclasses import dataclass
from enum import Enum, auto
from typing import List, Optional, Sequence, Tuple
import pathlib


# ----------------------------------------------------------------------
# 1. 基本枚举与数据结构
# ----------------------------------------------------------------------


class SystemType(Enum):
    """双透镜结构类型：胶合 / 分离式"""

    CEMENTED = auto()
    AIR_SPACED = auto()


@dataclass
class WavelengthSet:
    """
    工作谱线三元组（不再写死为 C-d-F）：

      lambda_ref   : 参考 / 主工作波长（原来的 d 线）
      lambda_short : 短波（原来的 F 线）
      lambda_long  : 长波（原来的 C 线）

    在可见光设计时你可以仍然用 (d, F, C)，
    在 UV / NIR / IR 设计时则换成你需要的三条线。
    """

    lambda_ref: float
    lambda_short: float
    lambda_long: float


@dataclass
class Glass:
    """
    玻璃库中的单块玻璃。

    注意命名：
    - n_catalog, V_catalog 是“玻璃目录给的参考 nd / Vd”，
      只是厂家在 C-d-F 下的标称值；
    - 真正用于当前设计波段的 Abbe 数由 abbe_number(glass, wavelengths)
      动态计算得到。
    """

    name: str
    n_catalog: float
    V_catalog: float
    formula: int
    coeffs: Tuple[float, ...] = ()  # 色散公式系数（CD 行），后续可用于 Sellmeier 等


@dataclass
class InputParams:
    """
    对应流程图顶部输入框中的参数。
    """

    f_prime: float
    D: float
    P0: float
    W0: float
    C0: float
    min_delta_V: float
    max_PE: float
    N: int
    system_type: SystemType
    wavelengths: WavelengthSet


@dataclass
class Geometry:
    """
    几何参数：曲率半径 & 厚度。
    对于分离式双透镜，你可以把 d1,d2 分别看成两个玻璃块的厚度；
    如果需要空气间隔，可再加 d_air。
    """

    R1: float
    R2: float
    R3: float
    d1: float
    d2: float
    # 可选：air gap 等
    geometric_ok: bool = True  # 几何可行性检查结果


@dataclass
class Invariants:
    """
    像差主参数 + 每一片的近轴光焦度。

    注意：
    - P, W, C 是三阶像差主参数（球差 / 彗差 / 轴向色差），
      在本简化 demo 中，P/W/C 先占位，由 compute_metrics 重新计算。
    - phi1, phi2 是两片透镜的近轴光焦度，真正决定几何结构。
    """

    P: float  # spherical aberration parameter (system-level), 先占位
    W: float  # coma parameter, 先占位
    C: float  # axial chromatic parameter, 通常由 C0 指定或后面计算
    phi1: float  # 第一片透镜的近轴光焦度
    phi2: float  # 第二片透镜的近轴光焦度


@dataclass
class AchromatState:
    """
    存储由 C = C0 色差条件得到的中间结果。
    目前只包含两片透镜的“相对焦度”，即:
        φ1 = k * phi1_rel
        φ2 = k * phi2_rel
    其中标度因子 k 由后续 P/W 方程决定。
    """

    phi1_rel: float
    phi2_rel: float


@dataclass
class Metrics:
    """
    像差评价：LCA, S1, W, PE 等。
    流程图中明确出现的是：LCA, S1, W, 预评估 PE。
    """

    LCA: float
    S1: float
    W: float
    PE: float


@dataclass
class DoubletResult:
    """
    单个候选双透镜的完整结果，用于 ranking。
    """

    glass1: Glass
    glass2: Glass
    geom: Geometry
    invariants: Invariants
    metrics: Metrics


# ----------------------------------------------------------------------
# 2. 玻璃库与简单工具函数
# ----------------------------------------------------------------------


def read_glass_catalog(path: Optional[str] = None) -> List[Glass]:
    """
    读取 Zemax 的 .AGF 文件。
    若不提供路径，则自动读取当前目录下的 SCHOTT.AGF。
    """
    if path is None:
        path = pathlib.Path(__file__).parent / "SCHOTT.AGF"
    else:
        path = pathlib.Path(path)

    glasses: List[Glass] = []
    current: Optional[Glass] = None

    with path.open("r", encoding="latin-1") as f:
        for raw_line in f:
            line = raw_line.strip()
            if not line or line.startswith("!"):
                continue

            # New Material
            if line.startswith("NM "):
                parts = line.split()
                if len(parts) < 6:
                    continue
                try:
                    name = parts[1]
                    formula = int(parts[2])
                    n_cat = float(parts[4])  # 玻璃库参考折射率（通常是 nd）
                    V_cat = float(parts[5])  # 玻璃库参考 Abbe 数（通常是 Vd）
                except Exception:
                    continue

                current = Glass(
                    name=name,
                    n_catalog=n_cat,
                    V_catalog=V_cat,
                    formula=formula,
                )
                glasses.append(current)

            # Coefficients (CD)
            elif line.startswith("CD ") and current is not None:
                try:
                    coeffs = tuple(float(x) for x in line.split()[1:])
                    current.coeffs = coeffs
                except Exception:
                    pass

    return glasses


# ----------------------------------------------------------------------
# 3. 色散 / 折射率相关工具（你后面可以逐步补充）
# ----------------------------------------------------------------------


def refractive_index(glass: Glass, wavelength: float) -> float:
    """
    给定玻璃与波长，返回折射率 n(λ)。

    当前版本：占位实现，简单返回玻璃库参考折射率 n_catalog。
    一旦你按照 Zemax/论文实现公式（Sellmeier / Herzberger 等），
    只需要在这里用 glass.formula + glass.coeffs 解析即可，其他代码无需改动。

    TODO: 用玻璃的色散公式（formula, coeffs）实现真实的 n(λ)。
    """
    _ = wavelength  # 占位，避免未使用警告
    return glass.n_catalog


def abbe_number(glass: Glass, wls: WavelengthSet) -> float:
    """
    用当前设计使用的三条谱线，计算该玻璃的 Abbe 数：

        V = (n_ref - 1) / (n_short - n_long)

    默认使用 refractive_index(glass, λ) 计算 n(λ)。
    如果当前占位实现导致 n_short ≈ n_long，使得分母过小，
    则退化为使用玻璃库给的 V_catalog 作为近似值。
    """
    n_ref = refractive_index(glass, wls.lambda_ref)
    n_short = refractive_index(glass, wls.lambda_short)
    n_long = refractive_index(glass, wls.lambda_long)

    denom = n_short - n_long
    if abs(denom) < 1e-12:
        # 暂时退化为目录 Abbe 数；当实现真实色散后就不会走到这里
        return glass.V_catalog

    return (n_ref - 1.0) / denom


def compute_partial_dispersion(glass: Glass, wls: WavelengthSet) -> float:
    """
    计算玻璃在选定谱线下的部分色散（如果论文中有用到）。
    这里先给一个占位接口，方便后续扩展。
    """
    _ = (glass, wls)
    # TODO: 按照论文中定义的部分色散公式实现
    return 0.0


# ----------------------------------------------------------------------
# 4. C = C0 色差方程 (Achromat condition)
# ----------------------------------------------------------------------


def solve_achromat_condition(
    g1: Glass, g2: Glass, params: InputParams
) -> Optional[AchromatState]:
    """
    解色差方程 C = C0（任意三波长版本）。

    按 ITMO/Nguyen 文章中的显式公式：
        φ1_rel = V1 (1 + V1 * C0) / (V1 - V2)
        φ2_rel = 1 - φ1_rel

    这里的 V1, V2 不再是“Vd”，而是使用 abbe_number(...)，
    根据 params.wavelengths 所指定的三条谱线计算得到。
    φ1_rel, φ2_rel 是“单位总功率”下的相对光焦度 (φ1 + φ2 = 1)，
    后续再通过 φ_total = 1/f' 统一缩放。
    """

    V1 = abbe_number(g1, params.wavelengths)
    V2 = abbe_number(g2, params.wavelengths)

    # 阿贝数差太小 → 本来就不太适合做消色差，还会导致公式不稳定
    if abs(V1 - V2) < 1e-9:
        return None

    invV1 = 1.0 / V1
    invV2 = 1.0 / V2

    num = params.C0 - invV2
    den = invV1 - invV2
    phi1_rel = num / den
    phi2_rel = 1.0 - phi1_rel
    return AchromatState(phi1_rel, phi2_rel)


def _scale_phi_from_fprime(
    achro_state: AchromatState, params: InputParams
) -> Optional[Tuple[float, float]]:
    """
    把相对光焦度缩放到真实光焦度（利用系统焦距 f'）：

        φ_total = 1/f'
        k       = φ_total / (φ1_rel + φ2_rel)
        φ1      = k φ1_rel
        φ2      = k φ2_rel
    """
    if abs(params.f_prime) < 1e-12:
        return None

    phi_total = 1.0 / params.f_prime
    denom = achro_state.phi1_rel + achro_state.phi2_rel
    if abs(denom) < 1e-12:
        return None

    k = phi_total / denom
    return k * achro_state.phi1_rel, k * achro_state.phi2_rel


# ----------------------------------------------------------------------
# 5. P, W 方程：胶合 & 分离式
# ----------------------------------------------------------------------


def solve_power_system_cemented(
    achro_state: AchromatState,
    params: InputParams,
) -> Optional[Invariants]:
    """
    胶合双透镜：在已知 C = C0（通常 C0 = 0）给出的 φ1:φ2 比例的前提下，
    利用系统焦距 f'（而不是 P0）来确定两片透镜的实际近轴光焦度 φ1, φ2。

    物理意义：
    - C = 0 决定两片玻璃分担光焦度的比例，使轴向色差抵消；
    - f' 决定总光焦度 φ_total = 1/f'；
    - P, W, C 作为像差主参数，将在后续 compute_metrics 阶段由曲率 R1,R2,R3 计算，
      而不是在这里被“伪装成光焦度”。
    """

    phi1_rel = achro_state.phi1_rel
    phi2_rel = achro_state.phi2_rel

    denom = phi1_rel + phi2_rel
    if abs(denom) < 1e-12:
        # 比例导致总焦度无法定义（几何上不可行的组合）
        return None

    # 系统总光焦度（目前采用高斯近轴近似：φ_total = 1 / f'）
    if abs(params.f_prime) < 1e-12:
        # 防御：不允许 f' = 0
        return None

    phi_total = 1.0 / params.f_prime
    k = phi_total / denom

    phi1 = k * phi1_rel
    phi2 = k * phi2_rel

    # 这里的 P, W, C 是三阶像差主参数，目前尚未根据 R1,R2,R3 计算，
    # 所以先用设计目标值占位（P0, W0, C0），真正的像差在 compute_metrics 中重新算。
    P_sph = params.P0  # 目标球差，可用于后续筛选
    W_coma = params.W0  # 目标彗差
    C_ax = params.C0  # 目标色差（通常为 0）

    invariants = Invariants(
        P=P_sph,
        W=W_coma,
        C=C_ax,
        phi1=phi1,
        phi2=phi2,
    )
    return invariants


def solve_power_system_air_spaced(
    achro_state: AchromatState,
    params: InputParams,
) -> Optional[Invariants]:
    """
    分离式双透镜：需要同时满足
        - 总光焦度 φ_total = 1/f'
        - 以及 W = W0（第三阶彗差条件）

    真正实现要解两个未知的方程组（通常涉及间隔和形状因子），
    这里先留空位，后续再逐步实现。
    """
    raise NotImplementedError(
        "Air-spaced doublet solver (P=W system) is not implemented yet."
    )


# ----------------------------------------------------------------------
# 6. 几何构造：由不变量 → R1, R2, R3, d1, d2
# ----------------------------------------------------------------------


def compute_geometry_from_invariants(
    invariants: Invariants,
    g1: Glass,
    g2: Glass,
    params: InputParams,
) -> Geometry:
    """
    根据两片透镜的近轴光焦度 phi1, phi2，构造一个物理合理的
    R1, R2, R3, d1, d2 组合（简化薄透镜近似版）。

    假设：
    - 第一片（g1）表面：R1（入射侧）、R2（胶合面）
    - 第二片（g2）表面：R2（胶合面）、R3（像侧）
    - 周围介质为空气 n=1
    - 使用薄透镜公式：
          phi1 ≈ (n1 - 1) * (1/R1 - 1/R2)
          phi2 ≈ (n2 - 1) * (1/R2 - 1/R3)

    这里 n1,n2 使用当前设计主波长 lambda_ref 下的折射率，
    以保持任意波段设计的一致性。
    """

    # 折射率用当前设计主波长
    n1 = refractive_index(g1, params.wavelengths.lambda_ref)
    n2 = refractive_index(g2, params.wavelengths.lambda_ref)

    phi1 = invariants.phi1
    phi2 = invariants.phi2

    # 1) 选择胶合面 R2：这里用一个简单的经验值，
    #    比如 R2 = -0.5 * f'，代表一个向物方弯曲的胶合面。
    #    将来你可以把它变成参数或由 shape factor 决定。
    R2 = -0.5 * params.f_prime
    if abs(R2) < 1e-6:
        R2 = -1.0  # 防御，避免除零

    inv_R2 = 1.0 / R2

    # 2) 根据薄透镜公式反解 R1, R3
    #    phi1 = (n1 - 1) * (1/R1 - 1/R2)
    #    => 1/R1 = phi1 / (n1 - 1) + 1/R2
    if abs(n1 - 1.0) < 1e-8:
        # 玻璃折射率异常（不应该等于空气），标记为不可行
        return Geometry(R1=0, R2=0, R3=0, d1=0, d2=0, geometric_ok=False)

    inv_R1 = phi1 / (n1 - 1.0) + inv_R2
    if abs(inv_R1) < 1e-12:
        # R1 → 无穷大（平面），这里先简单认为不可行
        return Geometry(R1=0, R2=0, R3=0, d1=0, d2=0, geometric_ok=False)
    R1 = 1.0 / inv_R1

    #    phi2 = (n2 - 1) * (1/R2 - 1/R3)
    #    => 1/R3 = 1/R2 - phi2 / (n2 - 1)
    if abs(n2 - 1.0) < 1e-8:
        return Geometry(R1=0, R2=0, R3=0, d1=0, d2=0, geometric_ok=False)

    inv_R3 = inv_R2 - phi2 / (n2 - 1.0)
    if abs(inv_R3) < 1e-12:
        return Geometry(R1=0, R2=0, R3=0, d1=0, d2=0, geometric_ok=False)
    R3 = 1.0 / inv_R3

    # 3) 厚度选择：先给一个和焦距成比例的合理值
    #    将来可以改成基于机械约束和口径 D 的函数。
    d1 = max(0.02 * params.f_prime, 1.0)  # 不小于 1 mm
    d2 = max(0.02 * params.f_prime, 1.0)

    geometric_ok = check_geometric_feasibility(R1, R2, R3, d1, d2, params)
    return Geometry(R1=R1, R2=R2, R3=R3, d1=d1, d2=d2, geometric_ok=geometric_ok)


def check_geometric_feasibility(
    R1: float,
    R2: float,
    R3: float,
    d1: float,
    d2: float,
    params: InputParams,
) -> bool:
    """
    几何和工艺约束检查（按 2022 文中的思路）：
      - 厚度必须为正；
      - 对每个表面要求 |R| > D/2，避免出现尖锐薄边。
    """

    if d1 <= 0 or d2 <= 0:
        return False

    half_aperture = params.D / 2.0
    for R in (R1, R2, R3):
        if abs(R) <= half_aperture:
            return False

    return True


# ----------------------------------------------------------------------
# 7. 像差与评价：Pi, Wi, Ci, LCA, S1, W, PE
# ----------------------------------------------------------------------


def compute_surface_invariants(
    geom: Geometry,
    g1: Glass,
    g2: Glass,
    params: InputParams,
) -> Tuple[Sequence[float], Sequence[float], Sequence[float]]:
    """
    对应流程图中的：
        Compute Pi, Wi, Ci

    这里实现一个“物理趋势正确但简化”的版本：

    - 每一面定义：
        φ_s = (n_out - n_in) / R_s     （薄面光焦度）
      其中：
        surface 1: air -> g1,   R1
        surface 2: g1  -> g2,   R2
        surface 3: g2  -> air,  R3

    - 球差权重（第三阶球差趋势）~ φ_s * y^4
    - 彗差权重（第三阶彗差趋势）~ φ_s * y^3 * h
      其中 y = D/2（边缘瞳高度），h 为像面场高（≈ f' * θ）

    注意：
    - 这不是完整的 Slusarev 公式，只是为后续替换精确 P,W,C 提供结构。
    - 真正论文里的系数可以直接替换 P_k, W_k 的表达式。
    """

    # 基本几何量
    y_pupil = params.D / 2.0  # 边缘瞳高度
    # 给一个默认场角（弧度），后续可以放到 InputParams 里作为参数
    theta_field = 0.05  # ~2.9°，只是一个 demo 数值
    h_image = params.f_prime * theta_field

    # 各面折射率：同样用当前设计主波长
    n_air = 1.0
    n1 = refractive_index(g1, params.wavelengths.lambda_ref)
    n2 = refractive_index(g2, params.wavelengths.lambda_ref)

    surfaces = [
        (n_air, n1, geom.R1),
        (n1, n2, geom.R2),
        (n2, n_air, geom.R3),
    ]

    P_surfaces: list[float] = []
    W_surfaces: list[float] = []
    C_surfaces: list[float] = []

    for n_in, n_out, R in surfaces:
        if abs(R) < 1e-9:
            # 曲率异常，给一个极大像差，便于后面筛掉
            P_surfaces.append(1e6)
            W_surfaces.append(1e6)
            C_surfaces.append(0.0)
            continue

        phi_s = (n_out - n_in) / R

        # 简化版：球差权重 ~ φ_s * y^4 / n_in^2
        base_n = max(abs(n_in), 1.0)
        P_k = phi_s * (y_pupil**4) / (base_n**2)

        # 简化版：彗差权重 ~ φ_s * y^3 * h / n_in^2
        W_k = phi_s * (y_pupil**3) * h_image / (base_n**2)

        # 色差权重 C_k：这里先设为 0（因为 C=0 设计），
        # 将来可以用 n(λ_short) - n(λ_long) 等公式替换。
        C_k = 0.0

        P_surfaces.append(P_k)
        W_surfaces.append(W_k)
        C_surfaces.append(C_k)

    return P_surfaces, W_surfaces, C_surfaces


def compute_LCA_S1_W_from_surfaces(
    P_surfaces: Sequence[float],
    W_surfaces: Sequence[float],
    C_surfaces: Sequence[float],
    geom: Geometry,
    g1: Glass,
    g2: Glass,
    params: InputParams,
) -> Tuple[float, float, float]:
    """
    按文中结构：
        LCA = Σ h_i C_i
        S1  = Σ h_i P_i
        W   = Σ h_i W_i   （这里用同样的权重，保持结构一致）

    这里用一个简单的近轴几何近似定义 h_i：
        h1 = 0
        h2 = d1
        h3 = d1 + d2
    """

    _ = (g1, g2, params)  # 目前未用到，先占位

    if len(P_surfaces) != 3:
        # 防御：当前只实现 3 个表面的情形（胶合双透镜）
        hs = [0.0] * len(P_surfaces)
    else:
        hs = [0.0, geom.d1, geom.d1 + geom.d2]

    LCA = sum(h * c for h, c in zip(hs, C_surfaces))
    S1 = sum(h * p for h, p in zip(hs, P_surfaces))
    W = sum(h * w for h, w in zip(hs, W_surfaces))

    return LCA, S1, W


def compute_PE_cemented(P: float, W: float, R2: float) -> float:
    """
    胶合双透镜的预评估参数 PE(P, W, R2)，按 2022 文中 (2) 的思路实现：

        - PE 随 |P|、|W| 增加而增大；
        - PE 随 |R2| 增大而减小（大半径 → 高阶像差小）；
        - 不显式使用 Q，由 1/|R2|^2 代替。

    这里采用一个简单的等价形式：
        PE = (|P|^3 * |W|) / |R2|^2

    系数和确切指数你可以在将来对照 Slusarev 原文后再精调。
    """

    if abs(R2) < 1e-9:
        return float("inf")

    return (abs(P) ** 3 * abs(W)) / (abs(R2) ** 2)


def compute_PE_air_spaced_from_P_surfaces(
    P_surfaces: Sequence[float],
) -> float:
    """
    分离式双透镜的 PE(P1,...,P4) 按 2022 文中公式 (3) 实现：

        PE = (1/4) * Σ |P_i|

    如果表面数不是 4，就按平均绝对值处理。
    """
    if not P_surfaces:
        return 0.0

    return sum(abs(p) for p in P_surfaces) / len(P_surfaces)


def compute_metrics(
    geom: Geometry,
    g1: Glass,
    g2: Glass,
    params: InputParams,
) -> Metrics:
    """
    封装以下步骤：
      - Compute Pi, Wi, Ci （每个表面）
      - LCA, S1, W （系统级）
      - PE （按 2022 文中的 (2) / (3) 区分胶合 / 分离）
    """
    P_surfaces, W_surfaces, C_surfaces = compute_surface_invariants(
        geom, g1, g2, params
    )
    LCA, S1, W = compute_LCA_S1_W_from_surfaces(
        P_surfaces, W_surfaces, C_surfaces, geom, g1, g2, params
    )

    if params.system_type == SystemType.CEMENTED:
        PE = compute_PE_cemented(P=S1, W=W, R2=geom.R2)
    else:
        PE = compute_PE_air_spaced_from_P_surfaces(P_surfaces)

    return Metrics(LCA=LCA, S1=S1, W=W, PE=PE)


# ----------------------------------------------------------------------
# 8. 单个玻璃对 (i,j) 的完整构建流程
# ----------------------------------------------------------------------


def dispersion_ok(g1: Glass, g2: Glass, params: InputParams) -> bool:
    """
    对应流程图菱形：
        |V[i] - V[j]| >= Min(ΔV) ?

    这里的 V 使用当前设计波段下的 abbe_number(...)，
    而不是目录中的 V_catalog。
    """
    V1 = abbe_number(g1, params.wavelengths)
    V2 = abbe_number(g2, params.wavelengths)
    return abs(V1 - V2) >= params.min_delta_V


def try_build_system_for_pair(
    g1: Glass,
    g2: Glass,
    params: InputParams,
) -> Optional[DoubletResult]:
    # 1) |V[i] - V[j]| >= Min(ΔV)?
    if not dispersion_ok(g1, g2, params):
        return None

    # 2) Solve C = C0
    achro_state = solve_achromat_condition(g1, g2, params)
    if achro_state is None:
        return None

    # 3) Solve "power system" (几何标度)，根据结构类型选择
    if params.system_type == SystemType.CEMENTED:
        invariants = solve_power_system_cemented(achro_state, params)
    else:
        invariants = solve_power_system_air_spaced(achro_state, params)

    if invariants is None:
        return None

    # 4) Compute geometry from invariants (用 phi1, phi2 推回 R1,R2,R3,d1,d2)
    geom = compute_geometry_from_invariants(invariants, g1, g2, params)
    if not geom.geometric_ok:
        return None

    # 5) Compute metrics (真正计算 P, W, C, LCA, S1, PE 等)
    metrics = compute_metrics(geom, g1, g2, params)
    if metrics.PE > params.max_PE:
        return None

    return DoubletResult(
        glass1=g1,
        glass2=g2,
        geom=geom,
        invariants=invariants,
        metrics=metrics,
    )


# ----------------------------------------------------------------------
# 9. ranking：维护前 N 个最优解 (cemented) / 直接保存 (air-spaced)
# ----------------------------------------------------------------------


def update_best_list(
    best_list: List[DoubletResult],
    candidate: DoubletResult,
    params: InputParams,
) -> None:
    """
    对应流程图下半部分：
      PE ≤ Max_PE? → Cemented objective? → size(doublet_res) < N? →
      Exists k: |W-W0| < |doublet_res[k].W-W0|? → 删除最差 → 保存当前系统

    为了可编辑性，这里尽量写成简单的线性逻辑。
    """

    # 非胶合：简单策略（论文中是直接保存或有轻微变化）
    if params.system_type != SystemType.CEMENTED:
        best_list.append(candidate)
        return

    # 胶合系统：Top-N by |W - W0|
    if len(best_list) < params.N:
        best_list.append(candidate)
        return

    worst_idx = None
    worst_dev = -1.0
    for idx, res in enumerate(best_list):
        dev = abs(res.metrics.W - params.W0)
        if dev > worst_dev:
            worst_dev = dev
            worst_idx = idx

    current_dev = abs(candidate.metrics.W - params.W0)
    if worst_idx is not None and current_dev < worst_dev:
        best_list[worst_idx] = candidate
    # 否则不做任何事（相当于流程图中的“Exists k...? → No”）


# ----------------------------------------------------------------------
# 10. 主循环：遍历玻璃对，综合所有步骤
# ----------------------------------------------------------------------


def synthesize_doublets(
    params: InputParams,
    glass_arr: Optional[Sequence[Glass]] = None,
) -> List[DoubletResult]:
    """
    对应整个流程图主循环：
      For i in glass_arr
        For j in glass_arr
          ...
      最后输出 doublet_res[] 列表。

    glass_arr:
      - 若为 None，则内部调用 read_glass_catalog()。
      - 否则使用传入的玻璃列表（方便你在 notebook 里测试子集）。
    """
    if glass_arr is None:
        glass_arr = read_glass_catalog()

    best_list: List[DoubletResult] = []

    for g1 in glass_arr:
        for g2 in glass_arr:
            candidate = try_build_system_for_pair(g1, g2, params)
            if candidate is None:
                continue
            update_best_list(best_list, candidate, params)

    # 最后整体排序（可选）：先按 |W-W0| 再按 PE
    best_list.sort(key=lambda r: (abs(r.metrics.W - params.W0), r.metrics.PE))
    return best_list

