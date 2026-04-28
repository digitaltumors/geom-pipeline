import shutil
import time
from pathlib import Path
from urllib.request import Request, urlopen

DOWNLOAD_TIMEOUT_SECONDS = 300.0
DOWNLOAD_RETRIES = 5
DOWNLOAD_BACKOFF_SECONDS = 5.0
DOWNLOAD_BUFFER_BYTES = 1024 * 1024 * 8


def download_to_path(source_url: str, output_path: Path) -> None:
	tmp_output = Path(f'{output_path}.part')
	request = Request(
		source_url,
		headers={
			'User-Agent': 'geom-pipeline-download/1.0',
		},
	)

	output_path.parent.mkdir(parents=True, exist_ok=True)

	for attempt in range(1, DOWNLOAD_RETRIES + 1):
		try:
			with (
				urlopen(request, timeout=DOWNLOAD_TIMEOUT_SECONDS) as response,
				tmp_output.open('wb') as handle,
			):
				shutil.copyfileobj(response, handle, length=DOWNLOAD_BUFFER_BYTES)
			tmp_output.replace(output_path)
		except Exception:
			tmp_output.unlink(missing_ok=True)
			if attempt == DOWNLOAD_RETRIES:
				raise
			sleep_seconds = DOWNLOAD_BACKOFF_SECONDS * attempt
			print(
				'[download_file] '
				f'attempt={attempt} failed; retrying in {sleep_seconds:.1f}s '
				f'url={source_url}',
				flush=True,
			)
			time.sleep(sleep_seconds)
		else:
			return


output_path = Path(snakemake.output[0])
source_url = snakemake.params.source_url
download_to_path(source_url, output_path)

print(
	'[download_file] '
	f'output={output_path} size_bytes={output_path.stat().st_size} '
	f'url={source_url} dataset_url={snakemake.params.get("dataset_url", "")}',
	flush=True,
)
