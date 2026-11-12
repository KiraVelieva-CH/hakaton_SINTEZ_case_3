from concurrent.futures import ThreadPoolExecutor, as_completed

def run_parallel(analyzers, context, workers=4):
    results, errors = {}, {}
    with ThreadPoolExecutor(max_workers=max(1, int(workers))) as pool:
        futures = {pool.submit(a.analyze, context): a.name for a in analyzers}
        for future in as_completed(futures):
            name = futures[future]
            try: results[name] = future.result()
            except Exception as e: errors[name] = str(e)
    return results, errors
