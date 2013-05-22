#ifndef CMANAGER_H
#define CMANAGER_H

/*
 * Class CComm provides necessary information for communication of data between two subdomains.
 *
 */
template <typename T>
class CComm
{
private:
	int _dstID;
	CVector<3,int> _send_size;
	CVector<3,int> _recv_size;
	CVector<3,int> _send_origin;
	CVector<3,int> _recv_origin;
	CVector<3,int> _normal;

public:
	CComm(int dstID):
	_dstID(dstID){

	}
	~CComm(){

	}
	CVector<3, int> getNormal() const {
		return _normal;
	}

	void setNormal(CVector<3, int> normal) {
		_normal = normal;
	}

	CVector<3, int> getRecvOrigin() const {
		return _recv_origin;
	}

	void setRecvOrigin(CVector<3, int> recvOrigin) {
		_recv_origin = recvOrigin;
	}

	CVector<3, int> getRecvSize() const {
		return _recv_size;
	}

	void setRecvSize(CVector<3, int> recvSize) {
		_recv_size = recvSize;
	}

	CVector<3, int> getSendOrigin() const {
		return _send_origin;
	}

	void setSendOrigin(CVector<3, int> sendOrigin) {
		_send_origin = sendOrigin;
	}

	CVector<3, int> getSendSize() const {
		return _send_size;
	}

	void setSendSize(CVector<3, int> sendSize) {
		_send_size = sendSize;
	}

	int getDstId() const {
		return _dstID;
	}

	void setDstId(int dstId) {
		_dstID = dstId;
	}

};
#endif
