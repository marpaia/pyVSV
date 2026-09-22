# Frontend

Next.js app. The API is `vsv.server` in `src/vsv/server/`.

From the repository root:

```bash
make dev
```

`make api` and `make frontend` start either process on its own.

`next dev` proxies `/api/*` to the API on port 8000. On Vercel, the root `vercel.json` sends `/api/*` to FastAPI and every other path to this app.
