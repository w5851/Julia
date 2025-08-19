"""
通用物理常数模块

此模块提供所有物理模型中使用的基本物理常数。
只包含真正通用的常数 - 模型特定参数在各自的模块中定义。
"""
module PhysicalConstants

export hc, Nc, MeV_to_fm_inv

# 通用物理常数
const hc = 197.33  # ℏc (MeV⋅fm) - 自然单位与SI单位间的转换因子
const Nc = 3       # QCD中的颜色数(通用规范理论常数)

# 单位转换函数
"""
    MeV_to_fm_inv(energy_MeV)

使用 ℏc = 197.33 MeV⋅fm 将能量从 MeV 转换为 fm⁻¹ 单位。
"""
MeV_to_fm_inv(energy_MeV) = energy_MeV / hc

end  # module PhysicalConstants
