using Pkg

println("开始安装项目依赖包到全局环境...")

# 根据 Project.toml 文件中的依赖项安装包
packages = [
    "ForwardDiff",
    "NLsolve", 
    "SpecialFunctions",
    "StaticArrays",
    "FiniteDifferences",
    "BenchmarkTools",
    "FastGaussQuadrature",
    "CSV",
    "PackageAnalyzer"  # 添加项目中需要的额外包
]

# 确保使用全局环境（不激活项目环境）
println("当前环境路径: ", Base.active_project())
if Base.active_project() !== nothing
    println("检测到项目环境已激活，正在切换到全局环境...")
    Pkg.activate()  # 激活全局环境
end
println("使用全局环境: ", Pkg.envdir())

# 安装每个包到全局环境
for pkg in packages
    println("正在安装包到全局环境: $pkg")
    try
        Pkg.add(pkg)
        println("✓ 成功安装到全局环境: $pkg")
    catch e
        println("✗ 安装失败: $pkg - $e")
    end
end

# 显示全局环境中的包状态
println("\n全局环境包状态:")
try
    Pkg.status()
catch e
    println("无法显示包状态: $e")
end

# 预编译包以提高后续加载速度
println("\n开始预编译包...")
try
    Pkg.precompile()
    println("✓ 预编译完成")
catch e
    println("✗ 预编译失败: $e")
end

println("\n" * "="^60)
println("安装完成！")
println("现在所有依赖包都已安装到全局环境中。")
println("可以直接运行脚本而无需激活项目环境。")
println("="^60)