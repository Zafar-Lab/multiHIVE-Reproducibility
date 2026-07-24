"""Plot runtime and memory scaling from the collected CSV results."""

from __future__ import annotations

from pathlib import Path

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import pandas as pd


SCRIPT_DIR = Path(__file__).resolve().parent
CSV_PATHS = {
	"multiHIVE": SCRIPT_DIR / "multiHIVE/memory_runtime_50epochs_results.csv",
	"multiVI": SCRIPT_DIR / "multiVI/memory_runtime_50epochs_results.csv",
}
OUTPUT_PATH = SCRIPT_DIR / "runtime_memory_plot.png"


def load_results(csv_path: Path, method: str) -> pd.DataFrame:
	dataframe = pd.read_csv(csv_path)
	
	dataframe = dataframe.copy()
	dataframe["Method"] = method
	dataframe["Num_Cells"] = pd.to_numeric(dataframe["Num_Cells"], errors="coerce")
	for column in ["Runtime_Seconds", "Peak_CPU_Memory_MB", "Peak_GPU_Memory_MB"]:
		if column not in dataframe.columns:
			dataframe[column] = pd.NA
		dataframe[column] = pd.to_numeric(dataframe[column], errors="coerce")
	dataframe["Runtime_Seconds"] = dataframe["Runtime_Seconds"]/50
	dataframe = dataframe.sort_values(["Method", "Dataset", "Num_Cells"])
	return dataframe


def plot_results(dataframe: pd.DataFrame, output_path: Path) -> None:
	metrics = [
		("Runtime_Seconds", "Runtime (seconds/epoch)"),
		("Peak_GPU_Memory_MB", "Peak GPU Memory (MB)"),
	]
	methods = list(dict.fromkeys(dataframe["Method"]))
	datasets = list(dict.fromkeys(dataframe["Dataset"]))
	colors = plt.get_cmap("tab10")
	linestyles = {
		"multiHIVE": "-",
		"multiVI": "--",
	}
	marker = {
		"multiHIVE": "o",
		"multiVI": "s",
	}

	fig, axes = plt.subplots(1, len(metrics), figsize=(13, 5), sharex=True)
	if len(metrics) == 1:
		axes = [axes]

	for axis, (metric, ylabel) in zip(axes, metrics):
		for dataset_index, dataset in enumerate(datasets):
			for method in methods:
				subset = dataframe[(dataframe["Method"] == method) & (dataframe["Dataset"] == dataset)]
				if subset.empty:
					continue
				axis.plot(
					subset["Num_Cells"],
					subset[metric],
					marker=marker.get(method, "o"),
					linewidth=2,
					color=colors(dataset_index),
					linestyle=linestyles.get(method, "-"),
					label=f"{dataset} ({method})",
				)
		axis.set_ylabel(ylabel)
		axis.grid(True, alpha=0.25)
		axis.set_xlabel("Number of cells")

	dataset_handles = [
		Line2D([0], [0], color=colors(index), linewidth=2, label=dataset)
		for index, dataset in enumerate(datasets)
	]
	method_handles = [
		Line2D([0], [0], color="black", linestyle=linestyles.get(method, "-"), linewidth=2, label=method)
		for method in methods
	]
	legend_dataset = fig.legend(
		handles=dataset_handles,
		title="Dataset",
		loc="upper left",
		bbox_to_anchor=(0.8, 0.7),
	)
	fig.add_artist(legend_dataset)
	fig.legend(
		handles=method_handles,
		title="Method",
		loc="lower left",
		bbox_to_anchor=(0.81, 0.3),
	)

	# fig.suptitle("Runtime and memory scaling by dataset", fontsize=16)
	fig.tight_layout(rect=(0, 0, 0.8, 0.97))
	fig.savefig(output_path, dpi=300, bbox_inches="tight")


def main() -> None:
	frames = [load_results(csv_path, method) for method, csv_path in CSV_PATHS.items()]
	dataframe = pd.concat(frames, ignore_index=True)
	plot_results(dataframe, OUTPUT_PATH)
	print(f"Saved plot to {OUTPUT_PATH}")


if __name__ == "__main__":
	main()
