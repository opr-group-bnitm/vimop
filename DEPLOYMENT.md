# ViMOP 本机部署记录

源码固定在 commit `1d72aa3221d84bb8881a7038198e309a7a8066d4`（ViMOP v1.0.6）。运行时使用 Nextflow 25.10.2；本机默认安装的 Nextflow 26.04.6 与上游配置不兼容，因此命令统一加 `NXF_VER=25.10.2`。

已完成：

- 官方数据库 `virus 2.4`、`contaminants 1.2`、`centrifuge 1.1` 下载、分片 SHA-256、合并归档 SHA-256、解压目录 SHA-256 全部通过。
- 数据库正式位置：`/data/databases/vimop/v1.0.6`，约 40 GB（2026-10-05 从 `/data/db/vimop` 迁入，三个子库 directory 级 SHA-256 复核通过）。
- 工作目录：`/data/work/vimop`；结果目录：`/gemini/pipelines/vimop/output`。
- Docker 已缓存全部 7 个镜像：`general:1.0.4`、`ingress:1.0.0`、`centrifuge:1.0.2`、`canu:1.0.1`、`medaka:1.0.3`、`report:1.0.0`、`structural_variants:1.0.0`（后三个经 mirror.gcr.io + aria2 分块下载 SHA-256 校验后 `docker load` 导入，缓存于 `/data/work/vimop/image_cache`，脚本 `pull_missing_images.py`）。

本地配置见 `deployment.config`。示例运行命令：

```bash
NXF_VER=25.10.2 nextflow run . -c deployment.config \
  --base_db /data/databases/vimop/v1.0.6 \
  --fastq /path/to/fastq \
  --out_dir /gemini/pipelines/vimop/output/sample
```

数据库下载的可断点校验脚本为 `download_official_db.py`，监测脚本为 `monitor.sh`。

## 验证运行（2026-10-05，已通过）

使用 `test_data/vimop-demo`（PBSIM3 模拟 LASV L/S + MS2 + 人源/棒杆菌背景，5000 条 ONT reads，位于 `/data/work/vimop/test_data/vimop-demo/lasv_simulated`）跑通完整流程，87/87 任务成功，耗时 10 分钟。三个预期基因组全部回收：

| 病毒 | 长度 | 深度 | 覆盖度 |
|---|---|---|---|
| LASV L（GU979513） | 7207 bp | 999x | 96.06% |
| LASV S（GU830839） | 3358 bp | 997x | 95.32% |
| MS2（LC710218） | 3605 bp | 798x | 93.31% |

结果在 `output/vimop-demo-validation/barcode01/`，含 `report_barcode01.html`。

已知注意事项：
- 本机 docker 为 containerd 存储，`docker load` OCI tar 后镜像以裸名（如 `medaka:1.0.3`）入库且该名无法直接 `docker tag`，须按镜像 ID 打标签（`pull_missing_images.py` 已处理）。
- `deployment.config` 的 `trace {}` 块会触发 Nextflow 25.10.2 TraceFileObserver NPE（不影响运行，trace.tsv 不更新）；介意可从配置中移除。
- 旧下载暂存 `/data/db/vimop/.staging-v1.0.6`（74G，含分卷+合并包+解包副本）在校验通过后已清理；`download_official_db.py` 与 `deployment.config` 的路径均已指向新位置。
