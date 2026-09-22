"use client";

import { useEffect, useState } from "react";

type Health = {
  status: string;
  service: string;
};

export default function Home() {
  const [health, setHealth] = useState<Health | null>(null);
  const [error, setError] = useState<string | null>(null);

  useEffect(() => {
    const controller = new AbortController();

    fetch("/api/health", { signal: controller.signal })
      .then(async (response) => {
        if (!response.ok) {
          throw new Error(`API responded ${response.status}`);
        }
        return (await response.json()) as Health;
      })
      .then((body) => {
        setHealth(body);
      })
      .catch((cause: unknown) => {
        if (cause instanceof DOMException && cause.name === "AbortError") {
          return;
        }
        setError(cause instanceof Error ? cause.message : "Could not reach the API");
      });

    return () => controller.abort();
  }, []);

  let status = "Checking the API…";
  if (health) {
    status = `${health.service} is ${health.status}`;
  } else if (error) {
    status = error;
  }

  return (
    <main className="mx-auto flex min-h-screen max-w-lg flex-col justify-center gap-3 px-6">
      <h1 className="text-3xl font-semibold tracking-tight">VSV</h1>
      <p className="text-neutral-600">{status}</p>
    </main>
  );
}
