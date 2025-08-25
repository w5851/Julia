using Pkg

println("开始安装项目依赖包...")

# 根据 Project.toml.backup 文件中的依赖项安装包
packages = [
    "ForwardDiff",
    "NLsolve",
    "SpecialFunctions", 
    "StaticArrays",
    "FiniteDifferences",
    "BenchmarkTools",
    "CSV"
]

# 激活当前项目环境
PROJECT_DIR = dirname(@__DIR__)  # 获取项目根目录
Pkg.activate(PROJECT_DIR)
println("激活项目环境: $PROJECT_DIR")

# 安装每个包
for pkg in packages
    println("正在安装包: $pkg")
    try
        Pkg.add(pkg)
        println("✓ 成功安装: $pkg")
    catch e
        println("✗ 安装失败: $pkg - $e")
    end
end

# 预编译包以提高后续加载速度
println("\n开始预编译包...")
try
    Pkg.precompile()
    println("✓ 预编译完成")
catch e
    println("✗ 预编译失败: $e")
end

println("\n安装完成！可以运行 `using Pkg; Pkg.status()` 查看已安装的包。")