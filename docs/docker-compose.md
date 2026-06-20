Yeah — interactive `docker run -it` plus a hand-typed `python` is the fragile part: if your SSH session drops the job can die, there's no restart, and the logs live only in your terminal. For a long single-machine job the clean setup is: run it **detached**, let the **wrapper be the container's command**, persist everything to a **mounted volume**, and manage the lifecycle declaratively. Docker Compose is the tidiest way to capture all of that.

One fix first — the wrapper currently calls `conda run -n reinvent4`, which buffers output and would make `docker logs -f` look dead. Let me switch it to activate the env directly so logs stream live and it works non-interactively.Now the Compose file that runs it detached with the right restart, stop-grace, and log settings:The lifecycle then becomes clean and disconnect-proof:

```bash
docker compose up -d --build     # build + start detached; survives SSH drop
docker compose logs -f           # live logs (streams now that conda run is gone)
docker compose stop              # graceful: SIGTERM -> checkpoint, waits up to 180s
docker compose up -d             # resume — wrapper picks up from the checkpoint
docker compose down              # stop + remove container (data is safe on the volume)
```

Why each piece matters: detached (`-d`) means the job is owned by the Docker daemon, not your terminal, so closing SSH doesn't touch it. `restart: on-failure` self-heals after a crash or OOM (the wrapper resumes from the checkpoint each time) but won't restart after a clean finish or a manual `stop` — Docker never auto-restarts a container you stopped by hand. `stop_grace_period: 180s` is the critical pairing with the SIGTERM checkpoint save: without it Docker SIGKILLs after 10s and you'd lose the in-flight chunk. The volume mount keeps configs, results, and checkpoints on the host so `docker compose down` and rebuilds don't lose anything.

A couple of operational notes. Monitoring is now just `docker compose logs -f` plus `docker stats reinvent4` for CPU/RAM and the usual `tail -f workspace/results/dockstream.log`. For TensorBoard, easiest is to run it on the host pointed at the mounted `workspace/.../tb_*` dir, or add a second tiny service in the same compose file mapping port 6006. And one caveat on `restart: on-failure`: if the job hits a *deterministic* failure (a config error, or it OOMs every time on the same chunk), it'll loop — resuming and re-failing — so glance at the logs after the first crash rather than assuming forward progress.

For right now, though, none of this applies to your validation pass — keep that one a plain foreground `python -m reinvent libinvent_rl.toml` with the small batch so you can watch it directly. Move to `docker compose up -d` once the loop is confirmed and you're launching the real multi-hundred-step run. If you tell me where your workspace lives on the host and whether this box is the CPU one or the new GPU machine, I'll fill in the exact volume path and the GPU block.
