# ViMOP CLAUDE.md

病毒基因組組裝流程（ONT 臨床/動物樣本，LASV/DENV/EBOV 等中小型 RNA 病毒最佳），OPR-group-BNITM v1.0.6。詳細部署記錄見同目錄 `DEPLOYMENT.md`；**新手學習/工具講解見 `docs/新手指南.md`（含彩色流程圖）**。

## 運行要點

- **必須 `NXF_VER=25.10.2`**：本機默認的 Nextflow 26.04.6 與上游配置不兼容。
- 標準命令：

```bash
NXF_VER=25.10.2 nextflow run /gemini/pipelines/vimop -c /gemini/pipelines/vimop/deployment.config \
  --fastq <fastq目錄或barcode上級目錄> \
  --out_dir /gemini/pipelines/vimop/output/<run名>
```

- work 目錄由 `deployment.config` 固定為 `/data/work/vimop`（NVMe SSD，符合全局規約）；結果輸出到項目 `output/`。
- 數據庫：`deployment.config` 中 `database_defaults.base = '/data/databases/vimop/v1.0.6'`（virus 2.4 / contaminants 1.2 / centrifuge 1.1，約 40 GB，2026-10-05 從 /data/db/vimop 遷入，校驅通過）。
- 7 個 Docker 鏡像（oprgroup/*）已全部導入本機；新增/重裝用 `pull_missing_images.py`（mirror.gcr.io + aria2 + SHA-256 校驗 + `docker load`）。

## 已知坑（本機實測）

- **docker load 後裸名無法 tag**：containerd 存儲下 OCI tar 導入的鏡像以裸名（`medaka:1.0.3`）顯示，`docker tag medaka:1.0.3 ...` 報 No such image；須先 `docker images --format '{{.Repository}}:{{.Tag}} {{.ID}}'` 拿 ID 再 tag。`pull_missing_images.py` 已內置此處理。
- **trace NPE**：`deployment.config` 的 `trace {}` 塊在 Nextflow 25.10.2 觸發 TraceFileObserver 空指針異常（刷屏但不影響運行，trace.tsv 不更新）。
- Docker Hub 直連極慢且最後層卡死；`mirror.gcr.io` 可正常拉 manifest 但 blob 建議走 aria2（腳本已實現）。
- 進程排查：nextflow 實際進程是 `java -jar nextflow-25.10.2-one.jar run ...`，`pgrep -f "nextflow run"` 匹配不到。

## 驗證狀態

2026-10-05 用 vimop-demo（LASV L/S + MS2 模擬數據，`/data/work/vimop/test_data/vimop-demo/lasv_simulated`）跑通完整流程，87/87 成功；三個預期基因組全部回收（覆蓋度 93-96%）。結果：`output/vimop-demo-validation/`。
