The neat part: TensorBoard only reads files, so a second container mounting the same workspace can serve the logs the running job is writing. Add a small service to your compose file:

```yaml
  tensorboard:
    image: htvs-pipeline:latest
    container_name: reinvent-tb
    command: ["/opt/conda/envs/reinvent4/bin/tensorboard",
              "--logdir", "/workspace/tb_libinvent_0",
              "--host", "0.0.0.0", "--port", "6006"]
    volumes:
      - ${PWD}:/workspace:ro
    ports:
      - "6006:6006"
    restart: unless-stopped
```

Then bring up just that service (it won't touch the running RL job):

```bash
docker compose up -d tensorboard
docker compose logs -f tensorboard      # confirm "TensorBoard ... at http://0.0.0.0:6006"
```

Or, if you'd rather not edit compose, a one-off does the same thing:

```bash
docker run --rm -d --name reinvent-tb -p 6006:6006 -v "$PWD:/workspace:ro" \
  htvs-pipeline:latest \
  /opt/conda/envs/reinvent4/bin/tensorboard \
  --logdir /workspace/tb_libinvent_0 --host 0.0.0.0 --port 6006
```

Three details that matter for the container version:

`--host 0.0.0.0` is mandatory inside a container. TensorBoard defaults to binding `localhost`, which means *localhost inside the container* — unreachable from outside even with the port published. Binding `0.0.0.0` is what makes `-p 6006:6006` actually work.

I mounted it `:ro` (read-only) — TensorBoard only reads, so there's no reason to give it write access to your run data.

Reaching the UI is the same as before: local machine → `http://localhost:6006`; remote/lab box → SSH tunnel from your laptop (`ssh -L 6006:localhost:6006 you@thehost`) then open `http://localhost:6006`.

And the same caveat as the host version once it's up: flip the Scalars time axis to RELATIVE or WALL so your one-step-per-chunk points spread over time instead of piling at step 0.

To stop it later: `docker compose stop tensorboard` (or `docker stop reinvent-tb` for the one-off). You can clean up the host venv too — `deactivate`
