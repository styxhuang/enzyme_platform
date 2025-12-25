from mxf_serum.ormp.entity import (BoxSelectParameter, Group, InputParameter,
                                   NodeConfigure, NodeHelpInfo, NodeIcon,
                                   NodeModel, NumberParameter,
                                   ParameterHelpInfo, SelectParameter,
                                   SwitchParameter)
from mxf_serum.ormp.postprocessing import Reflecting, SelectBoxEasy
from mxf_serum.ormp.relation import ChangeOptionsOn, ShowOn


class MacrosLigandsDockingParameter(NodeModel):
    __node__ = NodeConfigure(
        NodeHelpInfo(
            description='通过调用AutoDock Vina分子对接引擎执行多种分子对接计算任务，并可计算多种打分算法给出的打分值',
            publish_time='2022-10-23',
            references='https://vina.scripps.edu/',
            version='v1.3',
            software=dict(
                name='AutoDock Vina',
                version='1.2.3',
                license='Apache',
            ),
        ),
        '分子对接',
        'MOD_CADD_VINA_DOCKING_1_3.macros_ligand_docking',
        '利用分子对接引擎Autodock Vina执行分子对接',
        'LS-02-000-0003',
        component_icon=NodeIcon.READ_STRUCT_DATA,
        inport=2,
    )

    g1 = Group(
        '搜索区域',
        SelectParameter(
            'domain_type',
            '搜索区域查找方法',
            ParameterHelpInfo(
                description='设置查找搜索区域的方法',
            ),
            Reflecting({
                '自动查找': 1,
                '由对接盒子定义': 2,
            }),
            default='由对接盒子定义',
            options=['自动查找', '由对接盒子定义']
        ),
        SelectParameter(
            'auto_domain_method',
            '自动查找方法',
            ParameterHelpInfo(
                description='自动查找搜索区域的方法。FPocket利用基于泰森多边形的几何算法查找位点；AutoSite根据疏水性以及氢键供体受体等性质预测高亲和力位点。仅在“搜索区域查找方法”为自动查找时显示该项',
            ),
            ShowOn(depend_on='domain_type', depend_on_val='自动查找'),
            Reflecting({
                'FPocket': 'fpocket',
                'AutoSite': 'autosite',
            }),
            default='AutoSite',
            options=['FPocket', 'AutoSite'],
        ),
        NumberParameter(
            'auto_domain_num',
            '搜索区域个数',
            ParameterHelpInfo(
                description='自动查找搜索区域的个数最大值.仅在“搜索区域查找方法”为自动查找时显示该项',
            ),
            ShowOn(depend_on='domain_type', depend_on_val='自动查找'),
            default=1,
            vmin=1,
        ),
        NumberParameter(
            'auto_domain_lb',
            '搜索区域容积下限(A^3)',
            ParameterHelpInfo(
                description='',
            ),
            ShowOn(depend_on='domain_type', depend_on_val='自动查找'),
            default=0,
            precision=2,
            vmin=0,
        ),
        BoxSelectParameter(
            'box_select',
            '选择对接盒子',
            ParameterHelpInfo(
                description='通过组选择bounding box定义分子对接的搜索区域。仅在“搜索区域查找方法”为由对接盒子定义时显示该项',
            ),
            ShowOn(depend_on='domain_type', depend_on_val='由对接盒子定义'),
            SelectBoxEasy(),
        ),
    )

    g2 = Group(
        '对接参数',
        SelectParameter(
            'exhaustiveness',
            '搜索复杂度',
            ParameterHelpInfo(
                description='设置蒙特卡洛抽样次数，该参数越大，计算进行得更详细，可以得到更多的能量更好构象，但是也更耗时',
            ),
            default=8,
            options=[8, 16, 32, 64],
        ),
        SelectParameter(
            'scoring_func',
            '打分函数',
            ParameterHelpInfo(
                description='选择打分函数',
                recommendedValue='Vina',
            ),
            ChangeOptionsOn(
                depend_on='with_water',
                mapping={
                    True: ['Autodock4'],
                    False: [
                        'Autodock4',
                        'Vina',
                        'Vinardo',
                    ],
                }),
            Reflecting({
                'Autodock4': 'ad4',
                'Vina': 'vina',
                'Vinardo': 'vinardo',
            }),
            default='Vina',
            options=[
                'Autodock4',
                'Vina',
                'Vinardo',
            ],
        ),
        SwitchParameter(
            'flex_docking',
            '柔性对接',
            ParameterHelpInfo(
                description='是否启用柔性对接',
            ),
            default=False,
        ),
        InputParameter(
            'flex_res',
            '柔性侧链',
            ParameterHelpInfo(
                description='定义柔性对接中的柔性侧链, 仅在“柔性对接”为TRUE且“自动选择柔性侧链”为FALSE时显示该项',
            ),
            ShowOn(depend_on='flex_docking', depend_on_val=True)
        ),
    )

    g3 = Group(
        '构象输出控制',
        NumberParameter(
            'max_evals',
            '最大迭代次数',
            ParameterHelpInfo(
                description='设置进行分子对接搜索的最大迭代次数，后续会根据其他构象输出控制去筛选所获得的搜索结果。要注意的是，该值设置过低可能会导致所有搜索结果都被其他筛选条件过滤掉而导致无输出。该值为0则使用启发式搜索。',
                recommendedValue='0',
            ),
            default=0,
            vmin=0,
        ),
        NumberParameter(
            'min_rmsd',
            '最小RMSD差值',
            ParameterHelpInfo(
                description='构象间的最小RMSD差值（Angstrom）',
                recommendedValue='1.0',
            ),
            default=1.0,
            precision=1,
            step=0.1,
        ),
        NumberParameter(
            'energy_range',
            '能量范围',
            ParameterHelpInfo(
                description='与最佳构象的最大结合能差距（kcal/mol）',
                recommendedValue='3.0',
            ),
            default=3.0,
            precision=1,
            step=0.1,
        ),
        NumberParameter(
            'n_poses',
            '构象数量',
            ParameterHelpInfo(
                description='要输出的最大构象数量',
            ),
            default=9,
        ),
    )
