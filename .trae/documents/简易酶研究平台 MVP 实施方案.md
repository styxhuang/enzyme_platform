# 模块化容器架构 + K8s 可调度（MVP）实施方案

## 目标

* 按模块独立容器化（小分子、对接、动力学、分析、AlphaFold 接口）。

* 前后端与简单账号鉴权先行；打通“小分子结构预测模块”端到端。

* 共享宿主机文件夹挂载到各容器，K8s 下用 PV/PVC 统一存储，便于后台调度。

## 架构与容器

* 前端：`frontend`（React + Mol\* + Chart.js）。

* 后端：`api`（FastAPI + JWT + SQLite），统一作业调度与文件索引。

* 计算容器：

  * `mod-smol`（RDKit，小分子结构与性质）。

  * `mod-dock`（AutoDock Vina，对接）。

  * `mod-md`（GROMACS，动力学）。

  * `mod-ana`（分析指标）。

  * `mod-af2`（AlphaFold 接口占位）。

* 共享卷：开发用 Docker Compose `data:/data`；生产用 K8s `PersistentVolume`+`PersistentVolumeClaim` 挂载 `/data`。

## API 与作业模型

* 鉴权：`POST /auth/register`、`POST /auth/login`（JWT HttpOnly Cookie）、`POST /auth/logout`。

* 作业：`POST /jobs`（type+inputs）→ 后端路由到对应模块容器；`GET /jobs/{id}`、`GET /jobs`、`POST /jobs/{id}/retry`、`GET /files/{job_id}/{artifact}`。

* 模块内接口：统一 `POST /run|/predict`，约定写出到 `/data/<job_id>/outputs/` 并返回清单。

## 小分子结构预测（先实现）

* 输入：SMILES 列表/文件；参数：构象数、最小化方法（MMFF94/ UFF）。

* 流程：解析→ETKDG 构象→能量最小化→性质计算→导出（SDF/CSV/PNG 可选）。

* 输出指标：成功率、能量统计、性质分布；文件落地 `/data/<job_id>/outputs/`。

## 前端（现代简约、紧凑）

* 模块化导航：小分子、对接、动力学、分析、AlphaFold。

* 页面：登录/注册、仪表盘、各模块作业提交与结果列表、作业详情（下载工件、日志片段）。

* 视觉：浅色主题、紧凑控件、8pt 间距、响应式布局。

## K8s 调度与存储

* 每个模块以 Deployment/Service 暴露；后端以 ClusterIP 服务访问。

* 资源：`requests/limits` 明确，GPU 任务（未来 AF/MD）设置 `nodeSelector`/`resources`（如 `nvidia.com/gpu`）。

* 存储：PV/PVC（ReadWriteMany 优先，或 NFS/长卷），挂载 `/data` 到后端与计算容器。

## 实施步骤

1. 初始化仓库与目录：`frontend/`、`backend/`、`modules/`（smol/dock/md/analysis/af2）、`data/`、`docker-compose.yml`。
2. 后端：FastAPI 鉴权与作业模型；文件索引与下载接口；模块路由与占位实现。
3. 前端：登录与导航；作业提交表单（小分子）；列表与详情页骨架。
4. 小分子容器：RDKit 镜像与接口实现；打通端到端。
5. Compose 验证：本地多容器联调，确认共享卷读写。
6. K8s 清单（基础版）：Deployment/Service/Ingress（可选）、PV/PVC；在测试集群验证挂载与路由。
7. 后续：对接容器、分析占位、动力学容器占位、AlphaFold 接口页完善；事件推送与并发限流。

## 验收标准

* 登录后可提交小分子作业，查看状态并下载 SDF/CSV。

* 各模块页与接口占位可创建作业并落地规范化目录结构。

* K8s 部署验证：所有服务可发现；PV/PVC 挂载 `/data` 正常；作业链路打通。

