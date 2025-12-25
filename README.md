# enzyme_platform

## 启动

### 后端
- docker-compose up --build
- docker-compose build mod-md && docker-compose up -d --no-deps mod-md
- docker-compose build api && docker-compose up -d --no-deps api
- docker-compose build frontend && docker-compose up -d --no-deps frontend
- docker-compose build mod-dock && docker-compose up -d --no-deps mod-dock
- docker-compose build mod-mmpbsa && docker-compose up -d --no-deps mod-mmpbsa
- docker-compose build mod-analysis && docker-compose up -d --no-deps mod-analysis

mod-analysis
cd frontend && npm start dev

## 布局

按照以下描述进行优化

### 左边栏

顺序应为：
    - 小分子结构预测
    - 大分子结构预测 （原AlphaFold）
    - 分子对接 Docking
    - 分子动力学
    - 模拟结果分析
    - 任务列表

### 小分子预测

完成

### 大分子预测

左边：fasta序列，任务日志
右边：结构mol*

### 分子对接 Docking

- 受体：选择/上传，点击弹窗，选择结构或者上传
- 配体：输入SMILES，点击上传配体
- 结果列表：Vina的输出结果

### 分子动力学

1. 选择复合物，力场，水模型
2. 是否有小分子
    - 选择有小分子，选择/上传 小分子pdb，生成力场
    - 选择没有小分子，无动作
3. 右边的四块加上title，轨迹展示，结果表格，日志信息，结果图片
4. 展示图片需要考虑允许多张图片的切换
5. 现在刷新会自动跑到主页。但是动力学模拟本身是很慢的，所以要保证退出页面或者刷新都需要回到这个账户最后跑的任务，所以我想的是点击提交动力学之后，需要记录一个任务uuid，这样就会绑定一个任务。这个任务id需要缓存，每次进入到这个标签就加载这个任务的信息。
6. 现在只有在所有动力学都完成之后才会更新日志文件。但是我需要的是每个步骤都需要更新日志。日志的格式：正在进行EM，详细日志xxx.log，然后点击这个路径链接xxx.log，弹窗出现这个xxx.log。然后是EM已完成，开始NPT。
7. 需要在一个位置显示当前任务的uuid

### 模拟结果分析
1. 左边栏：
    1. 选择任务（可以选择动力学模拟过程的任务，通过任务的uuid）
    2. 选择性质分析（对最后的production的轨迹和能量进行分析，多选）
        - 通用性质分析（温度，压强，势能，动能， 点击选择后，出现后再打勾选择）
        - 结构性质分析
            - 选择组分（单选框，System， Protein， custom，如果选择了custom，显示字符框，用户填写index xx-xxx）
            - 选择性质分析（RMSD，RG，RMSF）
        - MMPBSA能量分析
            - 选择，MMPBSA的选项
    3. 性质分析按钮
    4. 日志框
2. 右边结果显示框：
    1. 右上角，加一个报告输出按钮
    2. 右边整体分成四部分：
        - 左上角：通用性质分析，显示图片，考虑多张图片切换
        - 左下角：结构性质分析，显示图片，考虑多张图片切换
        - 右上角：MMPBSA能量分析，显示图片，考虑多张图片切换
        - 右下角：结果表格：分成三个表1. 通用性质；2. 结构性质；3. MMPBSA分析

### MMPBSA能量计算

把mmpbsa单独变成一个模块，这个模块包含结合自由能计算，能量分解计算，丙氨酸扫描。其中先有一个分组设置

mmpbsa安装：
https://g-mmpbsa.readthedocs.io/en/latest/install.html
https://valdes-tresanco-ms.github.io/gmx_MMPBSA/dev/getting-started/

lt --port 8080 --subdomain biosim