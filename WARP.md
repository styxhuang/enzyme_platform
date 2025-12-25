# WARP.md

This file provides guidance to WARP (warp.dev) when working with code in this repository.

## Project Structure & Architecture

This is a multi-container scientific platform for enzyme research, consisting of a central backend, a frontend, and several specialized computational modules.

### Components

1.  **Frontend (`frontend/`)**
    *   **Type:** Node.js (Express) application serving static assets and UI.
    *   **Key Libs:** `molstar` (molecular visualization), `express`, `morgan`.
    *   **Entry Point:** `server.js`.

2.  **Backend (`backend/`)**
    *   **Type:** Python FastAPI application.
    *   **Role:** API Gateway, Authentication, Job Orchestration.
    *   **Database:** SQLite (`data/meta.sqlite`).
    *   **Entry Point:** `app.py`.
    *   **Key Logic:**
        *   Manages User Auth (JWT + SQLite).
        *   Dispatches jobs to module services via HTTP.
        *   Persists job state and outputs to disk and DB.

3.  **Modules (`modules/`)**
    *   Independent microservices built with FastAPI, each handling specific scientific tasks.
    *   **`smol`** (Port 8001): Small molecule processing (RDKit).
    *   **`af2`** (Port 8002): AlphaFold2 / Structure prediction.
    *   **`dock`** (Port 8003): Molecular docking (OpenBabel, Meeko).
    *   **`md`** (Port 8004): Molecular dynamics.
    *   **`analysis`** (Port 8005): Analysis tools.

### Data Flow

*   **Job Execution:** The backend receives a job request -> creates a record in SQLite -> calls the appropriate module API -> module processes data (reading/writing to shared `ENZYME_DATA_DIR`) -> backend updates job status.
*   **Storage:** Data is stored in a shared directory (default: `./data` locally, `/data` in containers), organized by Job ID or User ID.

## Development

### Prerequisites
*   Docker & Docker Compose
*   Node.js (for frontend local dev)
*   Python 3.11+ (for backend local dev)

### Common Commands

**Full Stack (Docker)**
```bash
# Build and start all services
docker-compose up --build
```

**Frontend**
```bash
cd frontend
npm install
# Run in development mode (auto-reload)
npm run dev
# Production start
npm start
```

**Backend**
```bash
cd backend
pip install -r requirements.txt
# Run locally (default port 8000)
uvicorn app:app --reload
```

### Testing
*   *Note: No automated test suite was detected in the repository structure.* 
*   When adding tests, prefer `pytest` for Python components and `jest` or `vitest` for the frontend.

### Coding Guidelines
*   **Python:** Follow standard PEP 8 guidelines. The backend and modules use FastAPI; ensure new endpoints use Pydantic models for validation.
*   **Database:** Schema changes in `backend/app.py` (`ensure_dirs` function) should be handled carefully as there is no formal migration tool (e.g., Alembic) currently set up; schema is applied on startup.
*   **Paths:** Always use `os.path.join` and reference `ENZYME_DATA_DIR` or `DATA_DIR` environment variables to ensure compatibility between local host and Docker environments.
